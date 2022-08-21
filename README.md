# Simple Lennard-Jones MD code optimized with MPI


## Usage
### makefile compiles the lattice generator "initCubicLattice.o" and the MD simulator "MD_mpi.o" -->
-   make
-   mpirun -np @{ranks} MD_mpi.o @{parameterfile}

### example:
-   generate parameter files specifying number of particles and density (in reduced units) with: ./initCubicLattice.o @{npart} @{density*}
-   run with: mpirun -np 4 MD_mpi.o /input/generated_grid.par

### output:
-   phys_ouput.dat > total, potential and kinetic energy for each timestep 
-   vtk_output >  vtk DataFile Version 4.0 for each participating rank 

## Versions
### gcc
- Apple clang version 13.1.6 (clang-1316.0.21.2.5) 
- Target: x86_64-apple-darwin21.5.0
- Thread model: posix
- InstalledDir: /Library/Developer/CommandLineTools/usr/bin
### mpi
- mpirun (Open MPI) 4.1.4
