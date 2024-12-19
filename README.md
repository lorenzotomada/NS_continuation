To run using MPI:

```bash
mpirun -np 4 ./main
```

References:
1) https://www.dealii.org/current/doxygen/deal.II/step_57.html (for Newton's method)
2) https://www.dealii.org/current/doxygen/deal.II/code_gallery_time_dependent_navier_stokes.html (for parallelization and time-dependence)
3) https://gitlab.com/lifex and in particular https://gitlab.com/lifex/lifex-cfd (for the SIMPLE preconditioner)
4) https://github.com/ICGonnella/SSFEM-Coanda-Effect/blob/main/source/coanda.cpp (for the 2D mesh)