/* -----------------------------------------------------------------------------
 TODO:
    1) Numerical discretization
      1.2) Ideally, theta-method (or at least K-N)
    2) Perform tests
      2.1) Find the steady, stable solution
      2.2) Without mesh refinement, symmetrical mesh
      2.3) Without mesh refinement, asymmetrical mesh
      2.4) With mesh refinement, symmetrical mesh
      2.5) With mesh refinement, symmetrical mesh, asymmetrical refinement (set refinement flag only e.g. if y > 7.5/2)
      2.6) Combine with initial guesses and see what changes (either, if possible, trying to obtain the other branch, or to get immediately the unstable one)
      2.7) At least hints for the 3D case
    3) Explain some details
      3.1) Why there is not a separate function for mesh-depending matrix and other components
      3.2) Why no trapezoidal rule
      3.3) Why need to keep first_iteration
* ------------------------------------------------------------------------------ */


#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/quadrature_point_data.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/parameter_handler.h>


#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/solver_selector.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>


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
#include <deal.II/numerics/nonlinear.h>


#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria.h>


#include <fstream>
#include <iostream>
#include <sstream>




// --------------------------------------------- NAMESPACE COANDA ---------------------------------------------



namespace coanda
{
  using namespace dealii;



// ------------------------------- TIME CLASS -------------------------------

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



// -------------------- TIME UTILITIES --------------------

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



// -------------------- DBCs --------------------

  template <int dim>
  class BoundaryValues : public Function<dim>
  {
  public:
    BoundaryValues() : Function<dim>(dim + 1) {}

    virtual double value(const Point<dim> &p, const unsigned int component) const override;
    virtual void vector_value(const Point<dim> &p, Vector<double> &values) const override;
  };



// -------------------- BOUNDARY VALUES VALUE FUNCTION --------------------
  
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
        double normalization_constant{0.5}; /* to leave the Reynolds number unchanged.
                                      The max inlet velocities doubles, so we need to rescale */
        value*=normalization_constant;
      }
      return value;
    }

    return 0;
  }



// -------------------- BOUNDARY VALUES VECTOR VALUE --------------------

  template <int dim>
  void BoundaryValues<dim>::vector_value(const Point<dim> &p, Vector<double> &values) const 
  {
    for (unsigned int c = 0; c < this->n_components; ++c) { values(c) = BoundaryValues<dim>::value(p, c); }
  }



// -------------------- SCHUR COMPLEMENT SPARSITY --------------------

  inline void schur_complement_sparsity(DynamicSparsityPattern &dst, const PETScWrappers::MPI::BlockSparseMatrix &src)
  {
    PETScWrappers::MPI::SparseMatrix tmp_mat;
    src.block(1, 0).mmult(tmp_mat, src.block(0, 1));
    //dst.reinit(tmp_mat.m(), tmp_mat.n());

    for (const auto &row : tmp_mat.locally_owned_range_indices())
    {
      for (auto entry = tmp_mat.begin(row); entry != tmp_mat.end(row); ++entry)
      {
        dst.add(entry->row(), entry->column());
      }
    }

    // Add entries for the (1, 1) block.
    for (const auto &row : src.block(1, 1).locally_owned_range_indices())
    {
      for (auto entry = src.block(1, 1).begin(row); entry != src.block(1, 1).end(row); ++entry)
      {
        dst.add(entry->row(), entry->column());
      }
    }
  }



// -------------------- SCHUR COMPLEMENT --------------------

  inline void schur_complement(PETScWrappers::MPI::SparseMatrix &dst, const PETScWrappers::MPI::BlockSparseMatrix &src, const PETScWrappers::MPI::Vector &diag_A00_inverse)
  {
    PETScWrappers::MPI::SparseMatrix tmp_mat;
    src.block(1, 0).mmult(tmp_mat, src.block(0, 1), diag_A00_inverse);
    dst = 0.0;
    dst.add(1.0, tmp_mat);
    dst.add(-1.0, src.block(1, 1));
  }





// -------------------- SIMPLE PRECONDITIONER --------------------

  class SIMPLE : public Subscriptor
  {
    public:
    SIMPLE(
          MPI_Comm mpi_communicator,
          TimerOutput &timer,
          const std::vector<IndexSet> &owned_partitioning,
          const std::vector<IndexSet> &relevant_partitioning,
          const PETScWrappers::MPI::BlockSparseMatrix &jacobian_matrix
    );


    void initialize_schur();
    void initialize();
    void assemble();

    void vmult(PETScWrappers::MPI::BlockVector &dst, const PETScWrappers::MPI::BlockVector &src) const
    {
      vmult_L(tmp, src);
      vmult_U(dst, tmp);
    }
  
    private:
    void vmult_L(PETScWrappers::MPI::BlockVector &dst, const PETScWrappers::MPI::BlockVector &src) const;
    void vmult_U(PETScWrappers::MPI::BlockVector &dst, const PETScWrappers::MPI::BlockVector &src) const;


    MPI_Comm mpi_communicator;
    TimerOutput &timer;

    const std::vector<IndexSet> &owned_partitioning;
    const std::vector<IndexSet> &relevant_partitioning;

    const SmartPointer<const PETScWrappers::MPI::BlockSparseMatrix> jacobian_matrix;
    PETScWrappers::MPI::SparseMatrix approximate_schur;

    mutable PETScWrappers::MPI::BlockVector tmp;
    PETScWrappers::MPI::Vector diag_F_inv;
  };



// --------------------- SIMPLE::CONSTRUCTOR ---------------------

  SIMPLE::SIMPLE(MPI_Comm mpi_communicator,
                TimerOutput &timer,
                const std::vector<IndexSet> &owned_partitioning,
                const std::vector<IndexSet> &relevant_partitioning,
                const PETScWrappers::MPI::BlockSparseMatrix &jacobian_matrix) : 
    mpi_communicator(mpi_communicator),
    timer(timer),
    owned_partitioning(owned_partitioning),
    relevant_partitioning(relevant_partitioning),
    jacobian_matrix(&jacobian_matrix)
  { 
    tmp.reinit(owned_partitioning, mpi_communicator);
  }



// --------------------- INITIALIZE SCHUR ---------------------

  void SIMPLE::initialize_schur()
  {
    DynamicSparsityPattern dsp_schur(relevant_partitioning[1]);
    schur_complement_sparsity(dsp_schur, (*jacobian_matrix));

    SparsityTools::distribute_sparsity_pattern(dsp_schur, owned_partitioning[1], mpi_communicator, relevant_partitioning[1]);
    approximate_schur.reinit(owned_partitioning[1], owned_partitioning[1], dsp_schur, mpi_communicator);
  }



// -------------------- INITIALIZE --------------------

  void SIMPLE::initialize()
  {
    diag_F_inv.reinit(owned_partitioning[0], mpi_communicator);

    tmp.reinit(owned_partitioning, mpi_communicator);
    //tmp = 0;
  }



// --------------------- ASSEMBLE ---------------------

  void SIMPLE::assemble()
  {
    // Initialize the sparsity pattern before the first call. This is done
    // here instead of in initialize() because the jacobian is not ready
    // before that.
    initialize_schur();
    initialize();

    const unsigned int start = diag_F_inv.local_range().first;
    const unsigned int end = diag_F_inv.local_range().second;

    for (unsigned int j = start; j < end; ++j)
    {
      AssertThrow((*jacobian_matrix).block(0, 0).diag_element(j) != 0.0, ExcDivideByZero());
      diag_F_inv(j) = 1.0 / (*jacobian_matrix).block(0, 0).diag_element(j);
    }

    diag_F_inv.compress(VectorOperation::insert);

    schur_complement(approximate_schur, (*jacobian_matrix), diag_F_inv);
  }



// --------------------- VMULT_L ---------------------

  void SIMPLE::vmult_L(PETScWrappers::MPI::BlockVector &dst, const PETScWrappers::MPI::BlockVector &src) const
  {
    TimerOutput::Scope timer_section(timer, "GMRES for preconditioner");

    SolverControl F_inv_control(20000, 1e-3 * src.block(0).l2_norm());
    PETScWrappers::SolverGMRES F_solver(F_inv_control, jacobian_matrix->get_mpi_communicator());

    PETScWrappers::PreconditionBoomerAMG F_inv_preconditioner(mpi_communicator, PETScWrappers::PreconditionBoomerAMG::AdditionalData());
    F_inv_preconditioner.initialize(jacobian_matrix->block(0, 0));

    F_solver.solve((*jacobian_matrix).block(0, 0), dst.block(0), src.block(0), F_inv_preconditioner); // d0 = F^{-1} * s0

    (*jacobian_matrix).block(1, 0).vmult(tmp.block(1), dst.block(0)); // t1 = (-B) * d0

    tmp.block(1).add(-1.0, src.block(1)); // t1 -= s1
    SolverControl Schur_inv_control(20000, 1e-3 * src.block(1).l2_norm());
    PETScWrappers::SolverGMRES Schur_solver(Schur_inv_control, mpi_communicator);

    PETScWrappers::PreconditionBoomerAMG Schur_inv_preconditioner(mpi_communicator, PETScWrappers::PreconditionBoomerAMG::AdditionalData());
    Schur_inv_preconditioner.initialize(approximate_schur);

    PETScWrappers::MPI::Vector solution(dst.block(1));

    solution = 0.0;
    Schur_solver.solve(approximate_schur, solution, tmp.block(1), Schur_inv_preconditioner); // d1 = (-Sigma^{-1}) * t1

    dst.block(1) = solution;

    //Schur_solver.solve(approximate_schur, dst.block(1), tmp.block(1), Schur_inv_preconditioner); // d1 = (-Sigma^{-1}) * t1
    // this throws an error: the initial guess can't be empty

    // In case of defective flow BC treatment the Navier-Stokes system matrix
    // has an extra block associated with the Lagrange multipliers (n_blocks=3).
    // Thus, in case of defective flow BCs we add an extra identity block.
    if (src.n_blocks() == 3) { dst.block(2) = src.block(2); }
  }



// --------------------- VMULT_U ---------------------

  void SIMPLE::vmult_U(PETScWrappers::MPI::BlockVector &dst, const PETScWrappers::MPI::BlockVector &src) const
  {
    dst.block(1) = src.block(1); // d1 = s1
    (*jacobian_matrix).block(0, 1).vmult(dst.block(0), src.block(1)); // d0 = B^T * s1
    dst.block(0).scale(diag_F_inv); // d0 = D^{-1} * d0
    dst.block(0).sadd(-1.0, src.block(0)); // d0 = s0 - d0

    if (src.n_blocks() == 3) { dst.block(2) = src.block(2); }
  }



// ------------------------------------ NS CLASS ----------------------------------- //

  template <int dim>
  class NS
  {
    public:
    NS(const bool adaptive_refinement,
      const bool use_continuation,
      const bool distort_mesh,
      const unsigned int fe_degree,
      const double stopping_criterion,
      const unsigned int n_glob_ref,
      const unsigned int jacobian_update_step,
      const double viscosity,
      const double viscosity_begin_continuation=0, // not necessary to pass this if continuation is not used
      const double continuation_step_size=0); // same

    void run();

    ~NS() { timer.print_summary(); }


    private:
     
    void make_grid();

    void setup_dofs();
    void setup_system();

    void assemble(const bool first_iteration, const bool assemble_jacobian, const bool steady_system=false); // assemble both jacobian and residual
    void assemble_system(const bool first_iteration, const bool steady_system=false);
    void assemble_rhs(const bool first_iteration, const bool steady_system=false);

    void solve(const bool first_iteration);

    void refine_mesh(const unsigned int min_grid_level, const unsigned int max_grid_level);

    void output_results(const unsigned int output_index, const bool output_steady_solution=false) const;

    void newton_iteration(PETScWrappers::MPI::BlockVector &dst,
                          const double tolerance,
                          const unsigned int max_n_line_searches,
                          const bool is_initial_step,
                          const bool steady_system=false);

    void compute_initial_guess();
    

    MPI_Comm mpi_communicator;
    const bool adaptive_refinement;
    const bool use_continuation;
    const bool distort_mesh; // done following https://www.dealii.org/current/doxygen/deal.II/step_49.html
    const unsigned int fe_degree; // added wrt step 57
    const double stopping_criterion; // same
    const unsigned int n_glob_ref;
    const unsigned int jacobian_update_step;
    double viscosity;
    const double viscosity_begin_continuation;
    const double continuation_step_size;

    parallel::distributed::Triangulation<dim> triangulation;
    FESystem<dim> fe;
    DoFHandler<dim> dof_handler;

    QGauss<dim> quad_formula;
    QGauss<dim-1> face_quad_formula;

    AffineConstraints<double> zero_constraints;
    AffineConstraints<double> nonzero_constraints;

    BlockSparsityPattern sparsity_pattern;

    PETScWrappers::MPI::BlockSparseMatrix jacobian_matrix;

    PETScWrappers::MPI::BlockVector present_solution;
    PETScWrappers::MPI::BlockVector old_solution;
    PETScWrappers::MPI::BlockVector steady_solution;

    PETScWrappers::MPI::BlockVector evaluation_points;
    PETScWrappers::MPI::BlockVector newton_update;
    PETScWrappers::MPI::BlockVector system_rhs;


    ConditionalOStream pcout;
    std::vector<types::global_dof_index> dofs_per_block;
    std::vector<IndexSet> owned_partitioning;
    std::vector<IndexSet> relevant_partitioning;
    IndexSet locally_relevant_dofs;

    Time time;
    mutable TimerOutput timer;

    std::shared_ptr<SIMPLE> preconditioner;

    std::vector<double> relative_distances;
    std::vector<double> residuals_at_first_iteration;
  };



// -------------------- CONSTRUCTOR --------------------

  template <int dim>
  NS<dim>::NS(const bool adaptive_refinement,
              const bool use_continuation,
              const bool distort_mesh,
              const unsigned int fe_degree,
              const double stopping_criterion,
              const unsigned int n_glob_ref,
              const unsigned int jacobian_update_step,
              const double viscosity,
              const double viscosity_begin_continuation,
              const double continuation_step_size
              )
    : mpi_communicator(MPI_COMM_WORLD),
      adaptive_refinement(adaptive_refinement),
      use_continuation(use_continuation),
      distort_mesh(distort_mesh),
      fe_degree(fe_degree),
      stopping_criterion(stopping_criterion),
      n_glob_ref(n_glob_ref),
      jacobian_update_step(jacobian_update_step),
      viscosity(viscosity),
      viscosity_begin_continuation(viscosity_begin_continuation),
      continuation_step_size(continuation_step_size),
      triangulation(mpi_communicator,
                    typename Triangulation<dim>::MeshSmoothing(Triangulation<dim>::smoothing_on_refinement |
                    Triangulation<dim>::smoothing_on_coarsening)),
      fe(FE_Q<dim>(fe_degree + 1), dim, FE_Q<dim>(fe_degree), 1),
      dof_handler(triangulation),
      quad_formula(fe_degree + 2),
      face_quad_formula(fe_degree + 2),
      pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0),
      time(20, 1e-2, 1e-1, 1e-1),
      timer(mpi_communicator, pcout, TimerOutput::never, TimerOutput::wall_times)  
  {
    make_grid();
  }



// -------------------- CREATE THE GRID --------------------

  template <int dim>
  void NS<dim>::make_grid()
  {
    Triangulation<dim> rectangle;   
   
    /* To create the mesh, we use subdivided_hyper_rectangle instead of hyper_rectangle because
    the latter yielded cells which were too streched. */

    if constexpr (dim == 2) // Checking with constexpr because it is known already at compile-time
    {
      std::vector<unsigned int> subdivisions{45, 10};
      GridGenerator::subdivided_hyper_rectangle(rectangle, subdivisions, Point<2>(0, 0), Point<2>(50, 7.5));
    }
    else
    {
      std::vector<unsigned int> subdivisions{50, 8, 8};
      GridGenerator::subdivided_hyper_rectangle(rectangle, subdivisions, Point<3>(0, 0, 0), Point<3>(50, 7.5, 7.5));
    }
    
    if (n_glob_ref > 0) { rectangle.refine_global(n_glob_ref); } // explain why it is here

    std::set<typename Triangulation<dim>::active_cell_iterator> cells_to_remove;
    bool inside_domain{true};
    for (const auto &cell : rectangle.active_cell_iterators())
    {
      inside_domain = true;
      for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
      {
        bool before_10_x_coord{cell->vertex(v)[0]<10};
        double v1{cell->vertex(v)[1]};
        bool first_check{(v1 > 5.0 && v1 < 7.5) || (v1 < 2.5)};
        if constexpr (dim == 2) 
        { 
          if (before_10_x_coord && first_check)
          {
            inside_domain = false; // if dim==2, it is a sufficient condition to be in the inlet
          }
        }
        else // otherwise we need to check also the z component
        {
          double v2{cell->vertex(v)[2]};
          bool second_check{(v2 > 5.0 && v2 < 7.5) || (v2 < 2.5)};
          if (before_10_x_coord && (first_check || second_check)) { inside_domain = false; }
        }
      }
      if (!inside_domain) { cells_to_remove.insert(cell); }
    }

    GridGenerator::create_triangulation_with_removed_cells(rectangle, cells_to_remove, triangulation);

    for (const auto &face : triangulation.active_face_iterators()) {
      if (face->at_boundary())
      {
        double face_center{face->center()[0]};
        if (std::fabs(face_center) < 1e-12)
        {
          face->set_boundary_id(1); // Inlet boundary
        }
        else
        {
          if (std::fabs(face_center - 50.0) < 1e-12)
          {
            face->set_boundary_id(2); // Outer boundary
          }
          else
          {
            face->set_boundary_id(3); // Wall boundary
          }
        }
      }
    }

    if (distort_mesh) { GridTools::distort_random(0.05, triangulation, true); }

    std::ofstream out("mesh.vtk");
    GridOut grid_out;
    grid_out.write_vtk(triangulation, out);
  }



// ------------------------------ SETUP DOFS ------------------------------

  template <int dim>
  void NS<dim>::setup_dofs()
  {
    dof_handler.distribute_dofs(fe);

    std::vector<unsigned int> block_component(dim + 1, 0);
    block_component[dim] = 1;
    DoFRenumbering::component_wise(dof_handler, block_component);
    dofs_per_block = DoFTools::count_dofs_per_fe_block(dof_handler, block_component);

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

    const FEValuesExtractors::Vector velocities(0);

    {
      nonzero_constraints.clear();
      zero_constraints.clear();

      nonzero_constraints.reinit(locally_relevant_dofs);
      zero_constraints.reinit(locally_relevant_dofs);

      DoFTools::make_hanging_node_constraints(dof_handler, nonzero_constraints);
      DoFTools::make_hanging_node_constraints(dof_handler, zero_constraints);

      // Apply Dirichlet boundary conditions on all Dirichlet boundaries except for the outlet.
      std::vector<unsigned int> dirichlet_bc_ids;
      dirichlet_bc_ids = std::vector<unsigned int>{1, 3};

      for (auto id : dirichlet_bc_ids)
      {
        VectorTools::interpolate_boundary_values(dof_handler,
                                                id,
                                                BoundaryValues<dim>(),
                                                nonzero_constraints,
                                                fe.component_mask(velocities)
                                                );
        VectorTools::interpolate_boundary_values(dof_handler,
                                                id,
                                                Functions::ZeroFunction<dim>(dim + 1),
                                                zero_constraints,
                                                fe.component_mask(velocities)
                                                );
      }

      nonzero_constraints.close();
      zero_constraints.close();
    }
  }

  
  
// --------------------- SETUP SYSTEM ---------------------
  
  template <int dim>
  void NS<dim>::setup_system()
  {
    preconditioner.reset();
    jacobian_matrix.clear();

    BlockDynamicSparsityPattern dsp(relevant_partitioning);
    DoFTools::make_sparsity_pattern(dof_handler, dsp, nonzero_constraints);
    sparsity_pattern.copy_from(dsp);

    SparsityTools::distribute_sparsity_pattern(dsp,
                                              dof_handler.locally_owned_dofs(),
                                              mpi_communicator,
                                              locally_relevant_dofs
                                              );
      
    jacobian_matrix.reinit(owned_partitioning, dsp, mpi_communicator);

    present_solution.reinit(owned_partitioning, relevant_partitioning, mpi_communicator);
    old_solution.reinit(owned_partitioning, relevant_partitioning, mpi_communicator);
    evaluation_points.reinit(owned_partitioning, relevant_partitioning, mpi_communicator);
    system_rhs.reinit(owned_partitioning, mpi_communicator); // non-ghosted
    newton_update.reinit(owned_partitioning, mpi_communicator); // non-ghosted
    steady_solution.reinit(owned_partitioning, relevant_partitioning, mpi_communicator);
  }



// -------------------- ASSEMBLE THE JACOBIAN --------------------

  template <int dim>
  void NS<dim>::assemble(const bool first_iteration, const bool assemble_jacobian, const bool steady_system)
  {
    TimerOutput::Scope timer_section(timer, "Assemble system");
    
    if (assemble_jacobian)
    {
      jacobian_matrix = 0;
    }

    system_rhs = 0;

    FEValues<dim> fe_values(fe, quad_formula, update_values | update_quadrature_points | update_JxW_values | update_gradients);
    
    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points = quad_formula.size();

    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar pressure(dim);

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     local_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    std::vector<Tensor<1, dim>> present_velocity_values(n_q_points);
    std::vector<Tensor<2, dim>> present_velocity_gradients(n_q_points);
    std::vector<Tensor<1, dim>> tmp(n_q_points); // present velocity - old solution velocity
    std::vector<double> present_velocity_divergences(n_q_points);
    std::vector<double> present_pressure_values(n_q_points);

    std::vector<double> div_phi_u(dofs_per_cell);
    std::vector<Tensor<1, dim>> phi_u(dofs_per_cell);
    std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell);
    std::vector<double> phi_p(dofs_per_cell);


    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
      {
        fe_values.reinit(cell);

        if (assemble_jacobian)
        {
          local_matrix = 0;
        }

        local_rhs = 0;

        fe_values[velocities].get_function_values(evaluation_points, present_velocity_values);
        fe_values[velocities].get_function_values(old_solution, tmp);
        fe_values[velocities].get_function_gradients(evaluation_points, present_velocity_gradients);
        fe_values[velocities].get_function_divergences(evaluation_points, present_velocity_divergences);
        fe_values[pressure].get_function_values(evaluation_points, present_pressure_values);

        for (size_t i = 0; i < tmp.size(); ++i) { tmp[i] -= present_velocity_values[i]; }

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
            if (assemble_jacobian)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
              {  
                double mass_matrix_contribution = ((phi_u[i] * phi_u[j]) / time.get_delta_t())
                                                   * fe_values.JxW(q);                //   From time stepping, 1/dt*M);
                local_matrix(i,j) += (
                      viscosity * scalar_product(grad_phi_u[i], grad_phi_u[j])        //   a(u, v) ---> mu*grad(u)*grad(v)  
                      - div_phi_u[i] * phi_p[j]                                       //   b(v, p) = -p*div(v)   
                      - phi_p[i] * div_phi_u[j]                                       //   b(u, q) = -q*div(u) 
                      + phi_u[i] * (present_velocity_gradients[q] * phi_u[j])         //   Linearization of the convective term (pt. 1) 
                      + phi_u[i] * (present_velocity_values[q] * grad_phi_u[j])       //   Linearization of the convective term (pt. 2) 
                      ) * fe_values.JxW(q);                                           //   JxW, integration weights   

                if (!steady_system) { local_matrix(i, j) += mass_matrix_contribution; }
              }
            }

            double present_velocity_divergence = trace(present_velocity_gradients[q]);

            double mass_matrix_contribution = ((phi_u[i]*tmp[q]) / time.get_delta_t())* fe_values.JxW(q);

            local_rhs(i) += (
                          -viscosity * scalar_product(grad_phi_u[i], present_velocity_gradients[q])
                          - phi_u[i] * (present_velocity_gradients[q] * present_velocity_values[q])
                          + div_phi_u[i] * present_pressure_values[q]
                          + phi_p[i] * present_velocity_divergence
                          ) * fe_values.JxW(q);

            if (!steady_system) { local_rhs(i) += mass_matrix_contribution; }
          }
        }


        cell->get_dof_indices(local_dof_indices);

        const AffineConstraints<double> &constraints_used = first_iteration ? nonzero_constraints : zero_constraints;

        if (assemble_jacobian)
        {
          constraints_used.distribute_local_to_global(local_matrix, local_rhs, local_dof_indices, jacobian_matrix, system_rhs);
        }
        else
        {
          constraints_used.distribute_local_to_global(local_rhs, local_dof_indices, system_rhs);
        }
      }
     
      if (assemble_jacobian)
      {
        jacobian_matrix.compress(VectorOperation::add);
        jacobian_matrix.block(1, 1) = 0;
      }

      system_rhs.compress(VectorOperation::add);
    }
  }



// -------------------- ASSEMBLE THE SYSTEM --------------------

  template <int dim>
  void NS<dim>::assemble_system(const bool first_iteration, const bool steady_system)
  {
    const bool assemble_jacobian{true};
    assemble(first_iteration, assemble_jacobian, steady_system);
  }
 


// -------------------- ASSEMBLE THE RESIDUAL --------------------

  template <int dim>
  void NS<dim>::assemble_rhs(const bool first_iteration, const bool steady_system)
  {
    const bool assemble_jacobian{false};
    assemble(first_iteration, assemble_jacobian, steady_system);
  }



// -------------------- SOLVE --------------------

  template <int dim>
  void NS<dim>::solve(const bool first_iteration)
  {
    const AffineConstraints<double> &constraints_used = first_iteration ? nonzero_constraints : zero_constraints;
 
    preconditioner.reset(new SIMPLE(mpi_communicator, timer, owned_partitioning, relevant_partitioning, jacobian_matrix));

    preconditioner -> assemble();

    SolverControl solver_control(/* jacobian_matrix.m() */ 50000, 1e-4 * system_rhs.l2_norm(), true);
    
    GrowingVectorMemory<PETScWrappers::MPI::BlockVector> vector_memory;
    SolverFGMRES<PETScWrappers::MPI::BlockVector> gmres(solver_control, vector_memory);

    gmres.solve(jacobian_matrix, newton_update, system_rhs, *preconditioner);
    // pcout << "FGMRES steps: " << solver_control.last_step() << std::endl;
 
    constraints_used.distribute(newton_update);
  }
 


// -------------------- NEWTON ITERATION --------------------

  template <int dim>
  void NS<dim>::newton_iteration(PETScWrappers::MPI::BlockVector &dst,
                                const double tolerance,
                                const unsigned int max_n_line_searches,
                                const bool is_initial_step,
                                const bool steady_system
                                )
  {
    bool first_iteration{is_initial_step};
    unsigned int line_search_n{0};
    double last_res{1.0};
    double current_res{1.0};
 
    while ((first_iteration || (current_res > tolerance)) && line_search_n < max_n_line_searches)
    {
      if (first_iteration)
      {
        setup_dofs();
        setup_system();

        evaluation_points = dst;

        assemble_system(first_iteration);
        solve(first_iteration);

        PETScWrappers::MPI::BlockVector tmp;
        tmp.reinit(owned_partitioning, mpi_communicator);
        tmp = newton_update; // initial condition is 0, no need to add anything (adding 0)
        
        nonzero_constraints.distribute(tmp);
        dst = tmp;
        // try to write  present_solution = newton_update;

        first_iteration = false;
        evaluation_points = dst;

        assemble_rhs(first_iteration);
        current_res = system_rhs.l2_norm();
        pcout << "  The residual of initial guess is " << current_res << std::endl;
        last_res = current_res;
      }
      else
      {
        // Choose the initial guess for Netwon's method: either the solution at the previous time step, or a linear combination 
        // between it and an asymmetrical one
        double guess_u_norm{steady_solution.block(0).l2_norm()};

        if (use_continuation && guess_u_norm!=0 && !steady_system && line_search_n==0)  // reinit initial guess
        {
          PETScWrappers::MPI::BlockVector tmp_initial_guess;
          tmp_initial_guess.reinit(owned_partitioning, mpi_communicator);

          tmp_initial_guess -= old_solution;
          tmp_initial_guess += steady_solution; // initial guess - old solution

          double dist_from_guess{tmp_initial_guess.l2_norm()};
          double alpha = dist_from_guess/guess_u_norm; // so that it is always in [0, 1]

          evaluation_points.reinit(owned_partitioning, relevant_partitioning, mpi_communicator);
          tmp_initial_guess.reinit(owned_partitioning, relevant_partitioning, mpi_communicator);
          tmp_initial_guess += old_solution; // if alpha = 0, || initial_guess - old solution || = 0, so that initial_guess = old solution
                                                 // and I should only weight it with the initial guess so alpha should multiply the previous solution
          tmp_initial_guess *= alpha;
          evaluation_points = steady_solution;  // (1-alpha)* initial_guess
          evaluation_points *= (1.0-alpha);
          evaluation_points += tmp_initial_guess; // alpha*old solution + (1-alpha)*initial_guess
        }

        else { evaluation_points = dst; }

        if (line_search_n%jacobian_update_step == 0 || line_search_n < 5) { assemble_system(first_iteration); /* We do not update the Jacobian at each iteration to reduce the cost */ }
        else { assemble_rhs(first_iteration); }
        solve(first_iteration);
 
        for (double alpha = 1.0; alpha > 1e-5; alpha *= 0.5)
        {
          evaluation_points.reinit(owned_partitioning, mpi_communicator);
          evaluation_points = dst;
          PETScWrappers::MPI::BlockVector tmp;
          tmp.reinit(owned_partitioning, mpi_communicator);
          tmp = newton_update;

          tmp *= alpha;
          evaluation_points += tmp;

          nonzero_constraints.distribute(evaluation_points);
          assemble_rhs(true); // we assemble taking into account the boundary conditions
          current_res = system_rhs.l2_norm();
          std::cout << "    alpha: " << std::setw(10) << alpha << std::setw(0) << "  residual: " << current_res << std::endl;
          if (current_res < last_res)
          {
            dst = evaluation_points;
            break;
          }
        }
                
        {
          std::cout << "  number of line searches: " << line_search_n << "  residual: " << current_res << std::endl;
          last_res = current_res;
          residuals_at_first_iteration.push_back(last_res);
        }
        
        ++line_search_n;
      }
    }
  }



// -------------------- OUTPUT RESULTS --------------------

  template <int dim>
  void NS<dim>::output_results(const unsigned int output_index, const bool output_steady_solution) const
  {
    TimerOutput::Scope timer_section(timer, "Output results");
    pcout << "Writing results..." << std::endl;

    std::vector<std::string> solution_names(dim, "velocity");
    solution_names.push_back("pressure");

    std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation(
                                                                              dim,
                                                                              DataComponentInterpretation::component_is_part_of_vector
                                                                              );
    data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);

    // vector to be output must be ghosted
    if (!output_steady_solution)
    { data_out.add_data_vector(present_solution, solution_names, DataOut<dim>::type_dof_data, data_component_interpretation); }

    else
    {
      PETScWrappers::MPI::BlockVector tmp;
      tmp.reinit(owned_partitioning, relevant_partitioning, mpi_communicator);
      tmp = steady_solution;
      data_out.add_data_vector(steady_solution, solution_names, DataOut<dim>::type_dof_data, data_component_interpretation); 
    }

    // Partition
    Vector<float> subdomain(triangulation.n_active_cells());
    
    for (unsigned int i = 0; i < subdomain.size(); ++i) { subdomain(i) = triangulation.locally_owned_subdomain(); }
    
    data_out.add_data_vector(subdomain, "subdomain");
    data_out.build_patches(fe_degree + 1);

    std::string basename;
    if (!output_steady_solution) { basename = "navierstokes" + Utilities::int_to_string(output_index, 6) + "-"; }
    else { basename = "steady_solution"; }

    std::string filename = basename + Utilities::int_to_string(triangulation.locally_owned_subdomain(), 4) + ".vtu";

    std::ofstream output(filename);
    data_out.write_vtu(output);

    static std::vector<std::pair<double, std::string>> times_and_names;
    
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      for (unsigned int i = 0; i < Utilities::MPI::n_mpi_processes(mpi_communicator); ++i)
      {
        times_and_names.push_back({time.current(), basename + Utilities::int_to_string(i, 4) + ".vtu"});
      }

      std::ofstream pvd_output(basename + ".pvd");
      DataOutBase::write_pvd_record(pvd_output, times_and_names);
    }
  }


 
// -------------------- REFINE MESH --------------------

  template <int dim>
  void NS<dim>::refine_mesh(const unsigned int min_grid_level, const unsigned int max_grid_level)
  {
    TimerOutput::Scope timer_section(timer, "Refine mesh");
    pcout << "Refining mesh..." << std::endl;

    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
    const FEValuesExtractors::Vector velocity(0);

    KellyErrorEstimator<dim>::estimate(dof_handler, face_quad_formula, {}, present_solution, estimated_error_per_cell, fe.component_mask(velocity));
    parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction(triangulation, estimated_error_per_cell, 0.6, 0.4);
    
    if (triangulation.n_levels() > max_grid_level)
    {
      for (auto cell = triangulation.begin_active(max_grid_level); cell != triangulation.end(); ++cell) { cell->clear_refine_flag(); }
    }

    for (auto cell = triangulation.begin_active(min_grid_level); cell != triangulation.end_active(min_grid_level); ++cell) { cell->clear_coarsen_flag(); }

    parallel::distributed::SolutionTransfer<dim, PETScWrappers::MPI::BlockVector> trans1(dof_handler);
    parallel::distributed::SolutionTransfer<dim, PETScWrappers::MPI::BlockVector> trans2(dof_handler);

    triangulation.prepare_coarsening_and_refinement();

    trans1.prepare_for_coarsening_and_refinement(present_solution);

    if (use_continuation)
    {
      trans2.prepare_for_coarsening_and_refinement(steady_solution); 
    }
    else { (void)trans2; }

    triangulation.execute_coarsening_and_refinement();

    // Reinitialize the system
    setup_dofs();
    setup_system();

    // Transfer solution
    // Need a non-ghosted vector for interpolation
    PETScWrappers::MPI::BlockVector tmp(newton_update);
    tmp = 0;
    trans1.interpolate(tmp);
    present_solution = tmp;

    if (use_continuation) 
    {
      tmp = 0;
      trans2.interpolate(tmp);
      steady_solution = tmp;
    }
  }



// -------------------- COMPUTE INITIAL GUESS --------------------

  template <int dim>
  void NS<dim>::compute_initial_guess()
  {
    bool is_initial_step = true;
    const double target_viscosity{viscosity};
 
    for (double mu = viscosity_begin_continuation; mu >= target_viscosity; mu -= continuation_step_size)
    {
      viscosity = mu;
      pcout << "Searching for initial guess with mu = " << mu << std::endl;

      const bool steady_system{true};
      newton_iteration(steady_solution, 1e-7, 50, is_initial_step, steady_system);
      
      is_initial_step = false;
      if (mu==target_viscosity) { break; }
    }

    viscosity = target_viscosity;
  }



// -------------------- RUN --------------------

  template <int dim>
  void NS<dim>::run()
  {
    pcout << "Running with PETSc on " << Utilities::MPI::n_mpi_processes(mpi_communicator) << " MPI rank(s)..." << std::endl;

    setup_dofs();
    setup_system();

    PETScWrappers::MPI::BlockVector tmp;
    tmp.reinit(owned_partitioning, mpi_communicator);

    nonzero_constraints.distribute(tmp);

    old_solution = tmp;

    // Time loop.

    bool should_stop{false};
    bool first_iteration{true};

    if (use_continuation)
    {
      compute_initial_guess();
      output_results(0, true); // save the steady solution
    }

    while ((time.end() - time.current() > 1e-12) && (!should_stop)) 
    {

      if (time.get_timestep() == 0) { output_results(0, false); }

      time.increment();
      std::cout.precision(6);
      std::cout.width(12);
      pcout << std::string(96, '*') << std::endl << "Time step = " << time.get_timestep() << ", at t = " << std::scientific << time.current() << std::endl;
      
      const bool steady_system{false};
      newton_iteration(present_solution, 1e-9, 75, first_iteration, steady_system);

      double norm_increment;

      {
        PETScWrappers::MPI::BlockVector solution_increment_in_time;
        solution_increment_in_time.reinit(owned_partitioning, mpi_communicator);
        solution_increment_in_time += present_solution;
        solution_increment_in_time -= old_solution;

        double old_u_norm = old_solution.block(0).l2_norm();
        double u_increment_norm = solution_increment_in_time.block(0).l2_norm();
        norm_increment = u_increment_norm/old_u_norm; // this could be inf, but no error is thrown
      }

      if (norm_increment < stopping_criterion) { should_stop = true; }
      pcout << " The relative distance between the two iterations is " << norm_increment << std::endl;
      relative_distances.push_back(norm_increment);

      old_solution.reinit(owned_partitioning, mpi_communicator);
      old_solution = present_solution;

      first_iteration = false;

      // Output
      if (time.time_to_output()) { output_results(time.get_timestep(), true); }

      if (time.time_to_refine())
      {
        if (adaptive_refinement) { refine_mesh(0, 2); }
      }
    }
  }
}





// ------------------------------- MAIN FUNCTION ------------------------------

int main(int argc, char *argv[])
{
  try
  {
    using namespace dealii;
    using namespace coanda;

    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
    
    const bool adaptive_refinement{true};
    const bool use_continuation{false};
    const bool distort_mesh{false};
    const unsigned int fe_degree{1};
    const double stopping_criterion{1e-10};
    const unsigned int n_glob_ref{1};
    const unsigned int jacobian_update_step{5};
    const double viscosity{0.7};
    const double viscosity_begin_continuation{0.72};
    const double continuation_step_size{1e-2};
    
    NS<2> bifurc_NS{adaptive_refinement,
                    use_continuation,
                    distort_mesh,
                    fe_degree,
                    stopping_criterion,
                    n_glob_ref,
                    jacobian_update_step,
                    viscosity,
                    viscosity_begin_continuation,
                    continuation_step_size};

    bifurc_NS.run();
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