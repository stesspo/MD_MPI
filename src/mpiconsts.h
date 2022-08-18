#ifndef MPICONSTS_H
#define MPICONSTS_H

#include <iostream>
#include <fstream>
#include <string.h>
#include <algorithm>
#include <mpi.h>
#include "consts.h"


typedef struct
{
    int x;
    int nx;
    MPI_Comm comm;
    // int neighbour_south;
    // int neighbour_north;
    int neighbour_east;
    int neighbour_west;
} Partition;

Partition createPartition(int mpi_rank, int mpi_size);
Domain createDomain(Partition p, Box boxdims);

#endif