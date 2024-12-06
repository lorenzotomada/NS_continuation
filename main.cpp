/* -----------------------------------------------------------------------------
 TODO: ideally in this order
    0) Change BCs of the initial condition ?
    1) Test for known values of mu (e.g. mu = 0.5)
    2) Netwon's method (+ preconditioner?)
    3) Relative distance between iterations as stopping criterion
    4) Various changes: pass variables to the constructor, maybe file for data as in step-35, ...
    5) Continuation algorithm!
    6) Find suitable initial guess for Newton's method (start from steady NS)
    7) What changes using mesh refinement?
 * ------------------------------------------------------------------------------ */

#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/quadrature_point_data.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparsity_tools.h>

#include <deal.II/lac/petsc_block_sparse_matrix.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria.h>

#include <fstream>
#include <iostream>
#include <sstream>

namespace coanda
{
  using namespace dealii;


  ////////// TIME CLASS //////////


  class Time
  {
  public:
    Time(const double time_end, const double delta_t, const double output_interval, const double refinement_interval)
      : timestep(0),
        time_current(0.0),
        time_end(time_end),
        delta_t(delta_t),
        output_interval(output_interval),
        refinement_interval(refinement_interval)
    {}

    double current() const { return time_current; }
    double end() const { return time_end; }
    double get_delta_t() const { return delta_t; }
    unsigned int get_timestep() const { return timestep; }

    bool time_to_output() const;
    bool time_to_refine() const;
    void increment();

  private:
    unsigned int timestep;
    double time_current;
    const double time_end;
    const double delta_t;
    const double output_interval;
    const double refinement_interval;
  };


  ///////// TIME UTILITIES //////////


  bool Time::time_to_output() const
  {
    unsigned int delta = static_cast<unsigned int>(output_interval / delta_t);
    return (timestep >= delta && timestep % delta == 0);
  }

  bool Time::time_to_refine() const
  {
    unsigned int delta = static_cast<unsigned int>(refinement_interval / delta_t);
    return (timestep >= delta && timestep % delta == 0);
  }

  void Time::increment()
  {
    time_current += delta_t;
    ++timestep;
  }


  ////////// DBCs //////////


  template <int dim>
  class BoundaryValues : public Function<dim>
  {
  public:
    BoundaryValues() : Function<dim>(dim + 1) {}

    virtual double value(const Point<dim> &p, const unsigned int component) const override;
    virtual void vector_value(const Point<dim> &p, Vector<double> &values) const override;
  };


  ////////// BOUNDARYVALUES VALUE FUNCTION //////////
  

  template <int dim>
  double BoundaryValues<dim>::value(const Point<dim> &p, const unsigned int component) const
  {
    Assert(component < this->n_components, ExcIndexRange(component, 0, this->n_components));
    double left_boundary = 0;
    if (component == 0 && std::abs(p[0] - left_boundary) < 1e-10)
    {
      double p1{p[1]};
      double value{20*(p1-2.5)*(5-p1)};
      if (dim == 3)
      {
        double p2{p[2]};
        value*=((p2-2.5)*(5-p2));
        double normalization_constant{0.5}; /* to leave the Reynolds number unchanged. The max inlet velocities doubles, so we need to rescale */
        value*=normalization_constant;
      }
      return value;
    }
    return 0;
  }




  template <int dim>
  void BoundaryValues<dim>::vector_value(const Point<dim> &p, Vector<double> &values) const 
  { for (unsigned int c = 0; c < this->n_components; ++c) { values(c) = BoundaryValues<dim>::value(p, c); } }

  
  class BlockSchurPreconditioner : public Subscriptor
  {
  public:
    BlockSchurPreconditioner(
      TimerOutput &timer,
      double gamma,
      double viscosity,
      double dt,
      const std::vector<IndexSet> &owned_partitioning,
      const PETScWrappers::MPI::BlockSparseMatrix &system,
      const PETScWrappers::MPI::BlockSparseMatrix &mass,
      PETScWrappers::MPI::BlockSparseMatrix &schur);

    void vmult(PETScWrappers::MPI::BlockVector &dst, const PETScWrappers::MPI::BlockVector &src) const;

  private:
    TimerOutput &timer;
    const double gamma;
    const double viscosity;
    const double dt;

    const SmartPointer<const PETScWrappers::MPI::BlockSparseMatrix> system_matrix;
    const SmartPointer<const PETScWrappers::MPI::BlockSparseMatrix> mass_matrix;
    const SmartPointer<PETScWrappers::MPI::BlockSparseMatrix> mass_schur;
  };

  BlockSchurPreconditioner::BlockSchurPreconditioner(
    TimerOutput &timer,
    double gamma,
    double viscosity,
    double dt,
    const std::vector<IndexSet> &owned_partitioning,
    const PETScWrappers::MPI::BlockSparseMatrix &system,
    const PETScWrappers::MPI::BlockSparseMatrix &mass,
    PETScWrappers::MPI::BlockSparseMatrix &schur)
    : timer(timer),
      gamma(gamma),
      viscosity(viscosity),
      dt(dt),
      system_matrix(&system),
      mass_matrix(&mass),
      mass_schur(&schur)
  {
    TimerOutput::Scope timer_section(timer, "CG for Sm");
    // The schur complemete of mass matrix is actually being computed here.
    PETScWrappers::MPI::BlockVector tmp1, tmp2;
    tmp1.reinit(owned_partitioning, mass_matrix->get_mpi_communicator());
    tmp2.reinit(owned_partitioning, mass_matrix->get_mpi_communicator());
    tmp1 = 1;
    tmp2 = 0;
    // Jacobi preconditioner of matrix A is by definition ${diag(A)}^{-1}$,
    // this is exactly what we want to compute.
    PETScWrappers::PreconditionJacobi jacobi(mass_matrix->block(0, 0));
    jacobi.vmult(tmp2.block(0), tmp1.block(0));
    system_matrix->block(1, 0).mmult(
      mass_schur->block(1, 1), system_matrix->block(0, 1), tmp2.block(0));
  }

  void BlockSchurPreconditioner::vmult(
    PETScWrappers::MPI::BlockVector &dst,
    const PETScWrappers::MPI::BlockVector &src) const
  {
    // Temporary vectors
    PETScWrappers::MPI::Vector utmp(src.block(0));
    PETScWrappers::MPI::Vector tmp(src.block(1));
    tmp = 0;
    // This block computes $u_1 = \tilde{S}^{-1} v_1$,
    // where CG solvers are used for $M_p^{-1}$ and $S_m^{-1}$.
    {
      TimerOutput::Scope timer_section(timer, "CG for Mp");
      SolverControl mp_control(src.block(1).size(), 1e-6 * src.block(1).l2_norm());
      PETScWrappers::SolverCG cg_mp(mp_control, mass_schur->get_mpi_communicator());
      // $-(\nu + \gamma)M_p^{-1}v_1$
      PETScWrappers::PreconditionBlockJacobi Mp_preconditioner;
      Mp_preconditioner.initialize(mass_matrix->block(1, 1));
      cg_mp.solve(
        mass_matrix->block(1, 1), tmp, src.block(1), Mp_preconditioner);
      tmp *= -(viscosity + gamma);
    }
    // $-\frac{1}{dt}S_m^{-1}v_1$
    {
      TimerOutput::Scope timer_section(timer, "CG for Sm");
      SolverControl sm_control(src.block(1).size(), 1e-6 * src.block(1).l2_norm());
      PETScWrappers::SolverCG cg_sm(sm_control, mass_schur->get_mpi_communicator());
      PETScWrappers::PreconditionNone Sm_preconditioner;
      Sm_preconditioner.initialize(mass_schur->block(1, 1));
      cg_sm.solve(mass_schur->block(1, 1), dst.block(1), src.block(1), Sm_preconditioner);
      dst.block(1) *= -1 / dt;
    }
    // Adding up these two, we get $\tilde{S}^{-1}v_1$.
    dst.block(1) += tmp;
    // Compute $v_0 - B^T\tilde{S}^{-1}v_1$ based on $u_1$.
    system_matrix->block(0, 1).vmult(utmp, dst.block(1));
    utmp *= -1.0;
    utmp += src.block(0);
    // Finally, compute the product of $\tilde{A}^{-1}$ and utmp
    // using another CG solver.
    {
      TimerOutput::Scope timer_section(timer, "CG for A");
      SolverControl a_control(src.block(0).size(), 1e-6 * src.block(0).l2_norm());
      PETScWrappers::SolverCG cg_a(a_control, mass_schur->get_mpi_communicator());
      PETScWrappers::PreconditionNone A_preconditioner;
      A_preconditioner.initialize(system_matrix->block(0, 0));
      cg_a.solve(
        system_matrix->block(0, 0), dst.block(0), utmp, A_preconditioner);
    }
  }


  ////////// NS CLASS //////////


  template <int dim>
  class NS
  {
  public:
    NS();
    void run();
    ~NS() { timer.print_summary(); }

  private:
    void make_grid();
    void setup_dofs();
    void make_constraints();
    void initialize_system();
    void assemble(bool use_nonzero_constraints, bool assemble_system);
    std::pair<unsigned int, double> solve(bool use_nonzero_constraints, bool assemble_system);
    void refine_mesh(const unsigned int, const unsigned int);
    void output_results(const unsigned int) const;
    int n_glob_ref;
    MPI_Comm mpi_communicator;
    parallel::distributed::Triangulation<dim> triangulation;
    double viscosity;
    double gamma;
    const unsigned int degree;

    std::vector<types::global_dof_index> dofs_per_block;

    FESystem<dim> fe;
    DoFHandler<dim> dof_handler;
    QGauss<dim> volume_quad_formula;
    QGauss<dim - 1> face_quad_formula;

    AffineConstraints<double> zero_constraints;
    AffineConstraints<double> nonzero_constraints;

    BlockSparsityPattern sparsity_pattern;
    PETScWrappers::MPI::BlockSparseMatrix system_matrix;
    PETScWrappers::MPI::BlockSparseMatrix mass_matrix;
    PETScWrappers::MPI::BlockSparseMatrix mass_schur;
    PETScWrappers::MPI::BlockVector present_solution;
    PETScWrappers::MPI::BlockVector solution_increment;
    PETScWrappers::MPI::BlockVector system_rhs;

    ConditionalOStream pcout;
    std::vector<IndexSet> owned_partitioning;
    std::vector<IndexSet> relevant_partitioning;
    IndexSet locally_relevant_dofs;
    std::shared_ptr<BlockSchurPreconditioner> preconditioner;

    Time time;
    mutable TimerOutput timer;
  };


  ////////// CONSTRUCTOR //////////


  template <int dim>
  NS<dim>::NS()
    : n_glob_ref(4),
      mpi_communicator(MPI_COMM_WORLD),
      triangulation(mpi_communicator, typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::smoothing_on_refinement | Triangulation<dim>::smoothing_on_coarsening)),
      viscosity(1.),
      gamma(0.1),
      degree(1),
      fe(FE_Q<dim>(degree + 1), dim, FE_Q<dim>(degree), 1),
      dof_handler(triangulation),
      volume_quad_formula(degree + 2),
      face_quad_formula(degree + 2),
      pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0),
      time(2*1e1, 1e-3, 1e-2, 1e-2),
      timer(mpi_communicator, pcout, TimerOutput::never, TimerOutput::wall_times)
  { make_grid(); }


  ////////// MAKE_GRID //////////


  template <int dim>
  void NS<dim>::make_grid(){
    Triangulation<dim> rectangle;
    if constexpr (dim == 2) { GridGenerator::hyper_rectangle(rectangle, Point<2>(0,0), Point<2>(50,7.5)); }
    else { GridGenerator::hyper_rectangle(rectangle, Point<3>(0,0,0), Point<3>(50, 7.5, 7.5)); }
    rectangle.refine_global(n_glob_ref);
    std::set<typename Triangulation<dim>::active_cell_iterator> cells_to_remove;
    bool inside_domain{true};
    for (const auto &cell : rectangle.active_cell_iterators())
    {
      inside_domain = true;
      for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
      {
        bool before_10_x_coord{cell->vertex(v)[0]<10};
        bool first_check{(cell->vertex(v)[1]>5 && cell->vertex(v)[1]<7.5) || (cell->vertex(v)[1]<2.5)};
        if constexpr (dim == 2) 
        { 
          if (before_10_x_coord && first_check) { inside_domain = false; /* if dim==2, it is a sufficient condition to be in the inlet */}
        }
        else
        {
          bool second_check{(cell->vertex(v)[2]>5 && cell->vertex(v)[2]<7.5) || (cell->vertex(v)[2]<2.5)};
          if (before_10_x_coord && (first_check || second_check)) { inside_domain = false; }
        }
      }
      if (!inside_domain) { cells_to_remove.insert(cell); }
    }
    GridGenerator::create_triangulation_with_removed_cells(rectangle, cells_to_remove, triangulation);
    for (const auto &face : triangulation.active_face_iterators()) {
      if (face->at_boundary())
      {
        if (std::fabs(face->center()[0]) < 1e-12) { face->set_boundary_id(1); /* Inlet boundary*/ }
        else
        {
          if (std::fabs(face->center()[0] - 50.0) < 1e-12) { face->set_boundary_id(2); /* Outer boundary */ }
          else { face->set_boundary_id(3); /* Wall boundary */ }
        }
      }
    }
  std::ofstream out("mesh_2d.vtk");
  GridOut grid_out;
  grid_out.write_vtk(triangulation, out);
}



  ////////// SETUP_DOFS //////////


  template <int dim>
  void NS<dim>::setup_dofs()
  {
    // The first step is to associate DoFs with a given mesh.
    dof_handler.distribute_dofs(fe);
    // We renumber the components to have all velocity DoFs come before
    // the pressure DoFs to be able to split the solution vector in two blocks
    // which are separately accessed in the block preconditioner.
    DoFRenumbering::Cuthill_McKee(dof_handler);
    std::vector<unsigned int> block_component(dim + 1, 0);
    block_component[dim] = 1;
    DoFRenumbering::component_wise(dof_handler, block_component);
    dofs_per_block = DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
    // Partitioning.
    unsigned int dof_u = dofs_per_block[0];
    unsigned int dof_p = dofs_per_block[1];
    owned_partitioning.resize(2);
    owned_partitioning[0] = dof_handler.locally_owned_dofs().get_view(0, dof_u);
    owned_partitioning[1] = dof_handler.locally_owned_dofs().get_view(dof_u, dof_u + dof_p);
    locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);
    relevant_partitioning.resize(2);
    relevant_partitioning[0] = locally_relevant_dofs.get_view(0, dof_u);
    relevant_partitioning[1] = locally_relevant_dofs.get_view(dof_u, dof_u + dof_p);
    pcout << "   Number of active fluid cells: " << triangulation.n_global_active_cells() << std::endl << "   Number of degrees of freedom: " << dof_handler.n_dofs() << " (" << dof_u << '+' << dof_p << ')' << std::endl;
  }

  // @sect4{InsIMEX::make_constraints}
  template <int dim>
  void NS<dim>::make_constraints()
  {
    // Because the equation is written in incremental form, two constraints
    // are needed: nonzero constraint and zero constraint.
    nonzero_constraints.clear();
    zero_constraints.clear();
    nonzero_constraints.reinit(locally_relevant_dofs);
    zero_constraints.reinit(locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, nonzero_constraints);
    DoFTools::make_hanging_node_constraints(dof_handler, zero_constraints);

    // Apply Dirichlet boundary conditions on all Dirichlet boundaries except for the outlet.
    std::vector<unsigned int> dirichlet_bc_ids;
    dirichlet_bc_ids = std::vector<unsigned int>{1, 3};

    const FEValuesExtractors::Vector velocities(0);
    for (auto id : dirichlet_bc_ids)
      {
        VectorTools::interpolate_boundary_values(dof_handler, id, BoundaryValues<dim>(), nonzero_constraints, fe.component_mask(velocities));
        VectorTools::interpolate_boundary_values(dof_handler, id, Functions::ZeroFunction<dim>(dim + 1), zero_constraints, fe.component_mask(velocities));
      }
    nonzero_constraints.close();
    zero_constraints.close();
  }


  ////////// INITIALIZE_SYSTEM //////////


  template <int dim>
  void NS<dim>::initialize_system()
  {
    preconditioner.reset();
    system_matrix.clear();
    mass_matrix.clear();
    mass_schur.clear();

    BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);
    DoFTools::make_sparsity_pattern(dof_handler, dsp, nonzero_constraints);
    sparsity_pattern.copy_from(dsp);
    SparsityTools::distribute_sparsity_pattern(dsp, dof_handler.locally_owned_dofs(), mpi_communicator, locally_relevant_dofs);

    system_matrix.reinit(owned_partitioning, dsp, mpi_communicator);
    mass_matrix.reinit(owned_partitioning, dsp, mpi_communicator);

    // Only the $(1, 1)$ block in the mass schur matrix is used.
    // Compute the sparsity pattern for mass schur in advance.
    // The only nonzero block has the same sparsity pattern as $BB^T$.
    BlockDynamicSparsityPattern schur_dsp(dofs_per_block, dofs_per_block);
    schur_dsp.block(1, 1).compute_mmult_pattern(sparsity_pattern.block(1, 0), sparsity_pattern.block(0, 1));
    mass_schur.reinit(owned_partitioning, schur_dsp, mpi_communicator);

    // present_solution is ghosted because it is used in the
    // output and mesh refinement functions.
    present_solution.reinit(owned_partitioning, relevant_partitioning, mpi_communicator);
    // solution_increment is non-ghosted because the linear solver needs
    // a completely distributed vector.
    solution_increment.reinit(owned_partitioning, mpi_communicator);
    // system_rhs is non-ghosted because it is only used in the linear
    // solver and residual evaluation.
    system_rhs.reinit(owned_partitioning, mpi_communicator);
  }

  // @sect4{InsIMEX::assemble}
  //
  // Assemble the system matrix, mass matrix, and the RHS.
  // It can be used to assemble the entire system or only the RHS.
  // An additional option is added to determine whether nonzero
  // constraints or zero constraints should be used.
  // Note that we only need to assemble the LHS for twice: once with the nonzero
  // constraint
  // and once for zero constraint. But we must assemble the RHS at every time
  // step.
  template <int dim>
  void NS<dim>::assemble(bool use_nonzero_constraints, bool assemble_system)
  {
    TimerOutput::Scope timer_section(timer, "Assemble system");

    if (assemble_system)
      {
        system_matrix = 0;
        mass_matrix = 0;
      }
    system_rhs = 0;

    FEValues<dim> fe_values(fe, volume_quad_formula, update_values | update_quadrature_points | update_JxW_values | update_gradients);
    FEFaceValues<dim> fe_face_values(fe, face_quad_formula, update_values | update_normal_vectors | update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points = volume_quad_formula.size();

    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar pressure(dim);

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> local_mass_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> local_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    std::vector<Tensor<1, dim>> current_velocity_values(n_q_points);
    std::vector<Tensor<2, dim>> current_velocity_gradients(n_q_points);
    std::vector<double> current_velocity_divergences(n_q_points);
    std::vector<double> current_pressure_values(n_q_points);

    std::vector<double> div_phi_u(dofs_per_cell);
    std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
    std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell);
    std::vector<double> phi_p(dofs_per_cell);

    for (auto cell = dof_handler.begin_active(); cell != dof_handler.end();
         ++cell)
      {
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);

            if (assemble_system)
              {
                local_matrix = 0;
                local_mass_matrix = 0;
              }
            local_rhs = 0;

            fe_values[velocities].get_function_values(present_solution, current_velocity_values);

            fe_values[velocities].get_function_gradients(present_solution, current_velocity_gradients);

            fe_values[velocities].get_function_divergences(present_solution, current_velocity_divergences);

            fe_values[pressure].get_function_values(present_solution, current_pressure_values);

            // Assemble the system matrix and mass matrix simultaneouly.
            // The mass matrix only uses the $(0, 0)$ and $(1, 1)$ blocks.
            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                for (unsigned int k = 0; k < dofs_per_cell; ++k)
                  {
                    div_phi_u[k] = fe_values[velocities].divergence(k, q);
                    grad_phi_u[k] = fe_values[velocities].gradient(k, q);
                    phi_u[k] = fe_values[velocities].value(k, q);
                    phi_p[k] = fe_values[pressure].value(k, q);
                  }

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  {
                    if (assemble_system)
                      {
                        for (unsigned int j = 0; j < dofs_per_cell; ++j)
                          {
                            local_matrix(i, j) +=
                              (viscosity *
                                 scalar_product(grad_phi_u[j], grad_phi_u[i]) -
                               div_phi_u[i] * phi_p[j] -
                               phi_p[i] * div_phi_u[j] +
                               gamma * div_phi_u[j] * div_phi_u[i] +
                               phi_u[i] * phi_u[j] / time.get_delta_t()) *
                              fe_values.JxW(q);
                            local_mass_matrix(i, j) +=
                              (phi_u[i] * phi_u[j] + phi_p[i] * phi_p[j]) *
                              fe_values.JxW(q);
                          }
                      }
                    local_rhs(i) -=
                      (viscosity * scalar_product(current_velocity_gradients[q],
                                                  grad_phi_u[i]) -
                       current_velocity_divergences[q] * phi_p[i] -
                       current_pressure_values[q] * div_phi_u[i] +
                       gamma * current_velocity_divergences[q] * div_phi_u[i] +
                       current_velocity_gradients[q] *
                         current_velocity_values[q] * phi_u[i]) *
                      fe_values.JxW(q);
                  }
              }

            cell->get_dof_indices(local_dof_indices);

            const AffineConstraints<double> &constraints_used = use_nonzero_constraints ? nonzero_constraints : zero_constraints;
            if (assemble_system)
              {
                constraints_used.distribute_local_to_global(local_matrix, local_rhs, local_dof_indices, system_matrix, system_rhs);
                constraints_used.distribute_local_to_global(local_mass_matrix, local_dof_indices, mass_matrix);
              }
            else { constraints_used.distribute_local_to_global(local_rhs, local_dof_indices, system_rhs); }
          }
      }

    if (assemble_system)
      {
        system_matrix.compress(VectorOperation::add);
        mass_matrix.compress(VectorOperation::add);
      }
    system_rhs.compress(VectorOperation::add);
  }

  // @sect4{InsIMEX::solve}
  // Solve the linear system using FGMRES solver with block preconditioner.
  // After solving the linear system, the same AffineConstraints object as used
  // in assembly must be used again, to set the constrained value.
  // The second argument is used to determine whether the block
  // preconditioner should be reset or not.
  template <int dim>
  std::pair<unsigned int, double>
  NS<dim>::solve(bool use_nonzero_constraints, bool assemble_system)
  {
    if (assemble_system)
      {
        preconditioner.reset(new BlockSchurPreconditioner(timer, gamma, viscosity, time.get_delta_t(), owned_partitioning,
                                                          system_matrix, mass_matrix, mass_schur));
      }

    SolverControl solver_control(
      system_matrix.m(), 1e-8 * system_rhs.l2_norm(), true);
    // Because PETScWrappers::SolverGMRES only accepts preconditioner
    // derived from PETScWrappers::PreconditionBase,
    // we use dealii SolverFGMRES.
    GrowingVectorMemory<PETScWrappers::MPI::BlockVector> vector_memory;
    SolverFGMRES<PETScWrappers::MPI::BlockVector> gmres(solver_control, vector_memory);

    // The solution vector must be non-ghosted
    gmres.solve(system_matrix, solution_increment, system_rhs, *preconditioner);

    const AffineConstraints<double> &constraints_used =
      use_nonzero_constraints ? nonzero_constraints : zero_constraints;
    constraints_used.distribute(solution_increment);

    return {solver_control.last_step(), solver_control.last_value()};
  }

  // @sect4{InsIMEX::run}
  template <int dim>
  void NS<dim>::run()
  {
    pcout << "Running with PETSc on " << Utilities::MPI::n_mpi_processes(mpi_communicator) << " MPI rank(s)..." << std::endl;
    //triangulation.refine_global(0);
    setup_dofs();
    make_constraints();
    initialize_system();

    // Time loop.
    bool refined = false;
    while (time.end() - time.current() > 1e-12)
      {
        if (time.get_timestep() == 0) { output_results(0); }
        time.increment();
        std::cout.precision(6);
        std::cout.width(12);
        pcout << std::string(96, '*') << std::endl << "Time step = " << time.get_timestep() << ", at t = " << std::scientific << time.current() << std::endl;
        // Resetting
        solution_increment = 0;
        // Only use nonzero constraints at the very first time step
        bool apply_nonzero_constraints = (time.get_timestep() == 1);
        // We have to assemble the LHS for the initial two time steps:
        // once using nonzero_constraints, once using zero_constraints,
        // as well as the steps immediately after mesh refinement.
        bool assemble_system = (time.get_timestep() < 3 || refined);
        refined = false;
        assemble(apply_nonzero_constraints, assemble_system);
        auto state = solve(apply_nonzero_constraints, assemble_system);
        // Note we have to use a non-ghosted vector to do the addition.
        PETScWrappers::MPI::BlockVector tmp;
        tmp.reinit(owned_partitioning, mpi_communicator);
        tmp = present_solution;
        tmp += solution_increment;
        present_solution = tmp;
        pcout << std::scientific << std::left << " GMRES_ITR = " << std::setw(3) << state.first << " GMRES_RES = " << state.second << std::endl;
        // Output
        if (time.time_to_output()) { output_results(time.get_timestep()); }
        if (time.time_to_refine())
          {
            refine_mesh(0, 4);
            refined = true;
          }
      }
  }

  // @sect4{InsIMEX::output_result}
  //
  template <int dim>
  void NS<dim>::output_results(const unsigned int output_index) const
  {
    TimerOutput::Scope timer_section(timer, "Output results");
    pcout << "Writing results..." << std::endl;
    std::vector<std::string> solution_names(dim, "velocity");
    solution_names.push_back("pressure");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    // vector to be output must be ghosted
    data_out.add_data_vector(present_solution, solution_names, DataOut<dim>::type_dof_data, data_component_interpretation);

    // Partition
    Vector<float> subdomain(triangulation.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
      { subdomain(i) = triangulation.locally_owned_subdomain(); }
    data_out.add_data_vector(subdomain, "subdomain");

    data_out.build_patches(degree + 1);

    std::string basename =
      "navierstokes" + Utilities::int_to_string(output_index, 6) + "-";

    std::string filename =
      basename + Utilities::int_to_string(triangulation.locally_owned_subdomain(), 4) + ".vtu";

    std::ofstream output(filename);
    data_out.write_vtu(output);

    static std::vector<std::pair<double, std::string>> times_and_names;
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        for (unsigned int i = 0;
             i < Utilities::MPI::n_mpi_processes(mpi_communicator);
             ++i)
          { times_and_names.push_back({time.current(), basename + Utilities::int_to_string(i, 4) + ".vtu"}); }
        std::ofstream pvd_output("navierstokes.pvd");
        DataOutBase::write_pvd_record(pvd_output, times_and_names);
      }
  }

  // @sect4{InsIMEX::refine_mesh}
  //
  template <int dim>
  void NS<dim>::refine_mesh(const unsigned int min_grid_level, const unsigned int max_grid_level)
  {
    TimerOutput::Scope timer_section(timer, "Refine mesh");
    pcout << "Refining mesh..." << std::endl;

    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
    const FEValuesExtractors::Vector velocity(0);
    KellyErrorEstimator<dim>::estimate(dof_handler, face_quad_formula, {}, present_solution, estimated_error_per_cell, fe.component_mask(velocity));
    parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction(
      triangulation, estimated_error_per_cell, 0.6, 0.4);
    if (triangulation.n_levels() > max_grid_level)
      {
        for (auto cell = triangulation.begin_active(max_grid_level); cell != triangulation.end(); ++cell)
          { cell->clear_refine_flag(); }
      }
    for (auto cell = triangulation.begin_active(min_grid_level); cell != triangulation.end_active(min_grid_level); ++cell)
      { cell->clear_coarsen_flag(); }

    // Prepare to transfer
    parallel::distributed::SolutionTransfer<dim, PETScWrappers::MPI::BlockVector>
      trans(dof_handler);

    triangulation.prepare_coarsening_and_refinement();

    trans.prepare_for_coarsening_and_refinement(present_solution);

    // Refine the mesh
    triangulation.execute_coarsening_and_refinement();

    // Reinitialize the system
    setup_dofs();
    make_constraints();
    initialize_system();

    // Transfer solution
    // Need a non-ghosted vector for interpolation
    PETScWrappers::MPI::BlockVector tmp(solution_increment);
    tmp = 0;
    trans.interpolate(tmp);
    present_solution = tmp;
  }
}





int main(int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      using namespace coanda;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      NS<2> flow{};
      flow.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl << "----------------------------------------------------" << std::endl;
      std::cerr << "Exception on processing: " << std::endl << exc.what() << std::endl << "Aborting!" << std::endl << "----------------------------------------------------" << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl << "----------------------------------------------------" << std::endl;
      std::cerr << "Unknown exception!" << std::endl << "Aborting!" << std::endl << "----------------------------------------------------" << std::endl;
      return 1;
    }
  return 0;
}