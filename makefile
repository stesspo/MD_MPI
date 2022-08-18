# by default select gcc
CFLAGS= -O3

SRC1 = ./src/main.cpp ./src/hostsources.cpp 
all: MD_mpi initGrid

MD_mpi:
	mpic++ $(SRC1) -o MD_mpi.o
initGrid:
	g++ ./src/initCubicLattice.cpp -o initCubicGrid.o 


clean:
	rm *.o vtk_output/*
