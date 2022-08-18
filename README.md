# Simple MD code optimized with MPI

Usage:

## makefile compiles the lattice generator "initCubicLattice.o" and the MD simulator "MD_mpi.o" -->
-   make
-   mpirun -np @{ranks} MD_mpi.o @{parameterfile}

## example:
-   generate parameter files specifying number of particles and density (in reduced units) with: ./initCubicLattice.o @{npart} @{density*}
-   run with: mpirun -np 4 MD_mpi.o /input/generated_grid.par


## The vtk_output folder stores files for visualisation in Paraview
