# Simple Lennard-Jones MD code optimized with MPI


## Usage
### makefile compiles the lattice generator "initCubicLattice.o" and the MD simulator "MD_mpi.o":
```bash
make
export OMP_NUM_THREADS = {ranks}
mpirun -np {ranks} MD_mpi.o {parameterfile}
```
### example:
-   generate parameter files specifying number of particles and density (in reduced units) with: ./initCubicLattice.o {npart} {density*} 
```bash
./initCubicLattice.o 1000 1
mpirun -np 4 MD_mpi.o /input/generated_grid.par
```
### output:
-   phys_ouput.dat > total, potential and kinetic energy for each timestep, can be plotted with
```bash
python plot.py
```
![Alt text](https://i.ibb.co/RHDtDhL/Mean-energy.png "Conservation of total Energy")

-   vtk_output >  vtk DataFile Version 4.0 for each participating rank 

Here a frame from ParaView-5.9.1 of two ranks (blue and red) communicating during interaction
![Alt text](https://i.ibb.co/tsxQxt4/blocksbigw.png "Blocks Big")

## Versions
### gcc
- Apple clang version 13.1.6 (clang-1316.0.21.2.5) 
- Target: x86_64-apple-darwin21.5.0
- Thread model: posix
### mpi
- mpirun (Open MPI) 4.1.4
