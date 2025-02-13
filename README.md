# Final project for the exam in Advanced Topics in Scientific Computing

The file `tutorial.md` constains a theoretical introduction to the material, as well as the commented code, and is structured in the style of a `deal.II` tutorial.

The file which is actually compiled when using CMake is, of course, `main.cpp`.

The folder `images` contains images needed for the tutorial, while the file `analysis.py` is just here for visualization purposes.

The file `separate_mesh_matrices.cpp` contains an implementation in which the matrices depending only on the mesh are implemented only after mesh refinement; however, for theoretical reasons explained in the tutorial, this is not a winning approach in this test case, and no tests using `separate_mesh_matrices.cpp` have been performed in addition to the one described in `tutorial.md`.

Instruction for compilation:
```bash
mkdir build && cd build

cmake -S ..
# wait for the process to finish

make
```

To run:
```bash
./main
```
or, using MPI:
```bash
mpirun -np 1 ./main
```
In the example above, `-np 1` is chosen to illustrate the use of the command.
However, for this particular problem, the moderate number of DoFs makes it inefficient to use multiple ranks, since overhead is unnecessarily introduced.

As a result, I recommend running the code in serial mode for testing purposes.
