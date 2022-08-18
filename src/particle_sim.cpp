#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string.h>
#include <cmath>
#include <mpi.h>
#include "hostsources.h"

using namespace std;

unsigned long get_time ()
{
    struct timeval tp;
    gettimeofday (&tp, NULL);
    return tp.tv_sec * 1000000 + tp.tv_usec;
}


int main(int argc, char* argv[])
{
  

  if(argc < 2){
		cout<<"Usage: mpirun -np {ranks} ./MD_mpi.o {name_file} \n \n";
    cout<<"- Example: \n ";
    cout<<"# generate a cubic lattice with: ./initCubicLattice.c {npart density*} \n \n";
    cout<<"# run with: mpirun -np 4  ./MD_mpi.o input/generated_grid.par \n \n";
		abort();
  }

  string file_par = argv[1];




// -------------------------------- MPI initialization
int provided;
MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
if (provided < MPI_THREAD_MULTIPLE) {
   cout<<"Error - MPI does not provide needed threading level"<<endl;
}
  // MPI_Init(&argc, &argv);
	int mpi_rank, mpi_size;

	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  int nparts_rank[4];


// Datatype used for exchanging particles information, it referes to:
// struct Field{
//   double x, y, z;
// }
  MPI_Datatype MPI_Field;
    int nblocks	=	1;
    int	blocklen[1] =	{3};	
    MPI_Datatype oldtypes[1]	=	{MPI_DOUBLE};	
    MPI_Aint displ[1]	=	{0};
    MPI_Type_create_struct(nblocks,	blocklen,	displ,	oldtypes,	&MPI_Field);	
    MPI_Type_commit(&MPI_Field);	
// --------------



// -------------- Parameters of experiment
  double timestep_length, time_end;
  double epsilon, sigma;
  double r_cut, r_skin;

  string* parameters = new string[24];

  int part_out_freq, vtk_out_freq;
  string input_file, part_out_name_base, vtk_out_name_base;
  int cl_workgroup_1dsize, cl_workgroup_3dsize_x, \
  cl_workgroup_3dsize_y, cl_workgroup_3dsize_z;
  float x_min, y_min, z_min, x_max, y_max, z_max;
  int x_n, y_n, z_n;
  int N;

loadparameters(file_par, parameters);

// read the parameter array
readparameters(parameters, input_file, time_end, timestep_length, epsilon, sigma, part_out_freq, part_out_name_base, vtk_out_freq, vtk_out_name_base, cl_workgroup_1dsize, cl_workgroup_3dsize_x, cl_workgroup_3dsize_y, cl_workgroup_3dsize_z, x_min, x_max, y_min, y_max, z_min, z_max, x_n, y_n, z_n, r_cut, r_skin);
delete[] parameters;

// Save box dimension
Box boxdim;
boxdim.x = x_max - x_min;
boxdim.y = y_max - y_min;
boxdim.z = z_max - z_min;
boxdim.volume = boxdim.x*boxdim.y*boxdim.z;
boxdim.r_min = r_skin + r_cut; //minimum cell length
// --------------------------------------------------------



// -------------- MPI: Generate 1D cartesian grid with east and west neighbours
 Partition p = createPartition(mpi_rank, mpi_size);
 Domain d = createDomain(p, boxdim);

  if(mpi_rank == 0) printf("Processor grid size %d\n", p.nx);
	printf("[Process %d]: Coordinate [%d]\n", mpi_rank, p.x);
  printf("[Process %d]: Neighbour E [%d], Neighbour W [%d]\n", mpi_rank, p.neighbour_east, p.neighbour_west);
  printf("[Process %d] Domain X: %f -> %f\n", mpi_rank, d.startx, d.endx);
  MPI_Barrier(MPI_COMM_WORLD);
// --------------------------------------------------------


// Read and store particles from input file, while send non-blocking infos about neighbours
// to involved processes. Here the initial memory allocation for neigh_p Field array is exact.
// Here is performed a first initialization to correctly allocate the initial memory to run the
// program, then the simulation can start.
Field *local_pos; // local positions
Field *local_vel; // local velocities
double *local_m; // local masses

get_input_particles(d, input_file, N, local_pos, local_vel, local_m, x_min, y_min, z_min);
Field *local_f = (Field*)malloc(d.npart*sizeof(Field));
printf("[Process %d]: number of local particles: %d\n", mpi_rank, d.npart);

// get size information from neighbouring ranks
h_count_neighbours(d, local_pos, boxdim, r_cut, r_skin);
get_neighbours_info(d, p, true);
d.memory_size = d.npart; //actual memory allocated into local arrays
// host buffers initialization


float d_add_mem = 1.3;  // additional memory added to arrays
// ----------------
Sizes buff_sizes;
buff_sizes.local_size = d.npart;
d.vacant = 0;

buff_sizes.size_neigh_send = int(d_add_mem*d.neigh_e_send) + int(d_add_mem*d.neigh_w_send);
buff_sizes.size_neigh_recv = int(d_add_mem*d.neigh_e_recv) + int(d_add_mem*d.neigh_w_recv);
buff_sizes.size_migrating_send = int(d_add_mem*d.migrating_e_send) + int(d_add_mem*d.migrating_w_send);
buff_sizes.size_migrating_recv = int(d_add_mem*d.migrating_e_recv) + int(d_add_mem*d.migrating_w_recv);
Field* h_neigh_sendBuffer = (Field*)malloc(buff_sizes.size_neigh_send*sizeof(Field));
Field* h_neigh_recvBuffer = (Field*)malloc(buff_sizes.size_neigh_recv*sizeof(Field));
Field* h_migrating_sendBuffer = (Field*)malloc(3*buff_sizes.size_migrating_send*sizeof(Field));
Field* h_migrating_recvBuffer = (Field*)malloc(3*buff_sizes.size_migrating_recv*sizeof(Field));
double* h_mass_sendBuffer = (double*)malloc(buff_sizes.size_migrating_send*sizeof(double));
double* h_mass_recvBuffer = (double*)malloc(buff_sizes.size_migrating_recv*sizeof(double));

int* mi; //index of migrating particles
mi = (int*) malloc (buff_sizes.size_migrating_send*sizeof(int));

MPI_Barrier(MPI_COMM_WORLD);

// send and recieve buffers
get_particles(d, p, local_pos, h_neigh_sendBuffer, h_neigh_recvBuffer, h_migrating_sendBuffer, h_migrating_recvBuffer, h_mass_sendBuffer, h_mass_recvBuffer, MPI_Field, true);
// -------------- 


// ---------------- iteration parameters:
int iters = (int) time_end/timestep_length + 1;
  int vtk = 0;
  int out = 0;
  int phys = 0;
  int phys_out_freq = 5;
// ----------------

// ---------------- Output initial conidtions
  ofstream phys_out;
if(p.x == 0) 
  phys_out.open("phys_output.dat");

 double en = 0;
 double v2 = 0;
 double global_en = 0;
 double global_k = 0;

  for(int i = 0; i < d.npart; i++){
    v2 += (pow(local_vel[i].x, 2) + pow(local_vel[i].y, 2) + pow(local_vel[i].z, 2))/(2*local_m[i]);
    local_f[i].x = 0; local_f[i].x = 0; local_f[i].z = 0;
    for(int j=0;j<i;j++){
      h_add2LJForce(local_f[i], local_f[j], local_pos[i], local_pos[j], en, epsilon, sigma, r_cut, r_skin, boxdim);
    }
    for(int k = 0;k < d.neigh_e_recv+d.neigh_w_recv;k++){
        h_addLJForce(local_f[i], local_pos[i], h_neigh_recvBuffer[k], en, epsilon, sigma, r_cut, r_skin, boxdim);
      }
  }
  MPI_Reduce(&en, &global_en, 1, MPI_DOUBLE, MPI_SUM, 0, p.comm);
  MPI_Reduce(&v2, &global_k, 1, MPI_DOUBLE, MPI_SUM, 0, p.comm);
// // E = Mechanical Energy, U = Potential Energy, K = Kinetic Energy, P = Total Momentum, T = Temperature
// if(p.x == 0){
//   phys_out<<"#time"<< setw(15) <<"E"<< setw(15) <<"U"<< setw(15) <<"K"<< setw(15) <<"P"<< setw(15) << "T"<<endl;
//   phys_out<< 0 << setw(10) << double((global_en + global_k)/N) << setw(15) << double(global_en/N) << setw(15) << double(global_k/N) <<endl;
// }
// ----------------


// ------------------ Migration flags
  bool flag_migration = false;
  int migration_freq = 1;

// ------------------------------------------ Start simulation
  if(p.x == 0){
    cout<<"Start Simulation..."<<endl;
  }

    MPI_Barrier(MPI_COMM_WORLD);
    double nTimeStop = 0;
    double nTimeLost = 0;
    double nTimeStart = get_time ();
// -----------------------------------------------------------
  h_LJForces(local_f, local_pos, h_neigh_recvBuffer, en, sigma, epsilon, r_cut, r_skin, boxdim, d.npart, d.neigh_e_recv+d.neigh_w_recv); 
// --------------------------- Start Loop


  for(int t = 0; t < iters; t++){
    en = 0;
    v2 = 0;
    if(t % 50 == 0 && p.x == 0)
      printf("iter: %i \n", t);
    
  
    if( t % migration_freq == 0)
      flag_migration = true;
    else
      flag_migration = false; 
    
    

    h_mv_step(local_pos, local_vel, local_f, local_m, timestep_length, boxdim, d.npart);
    h_count_neighbours(d, local_pos, boxdim, r_cut, r_skin);
    get_neighbours_info(d, p, flag_migration);

    h_check_dims(d, p, local_pos, local_vel, local_f, local_m, mi, h_neigh_sendBuffer, h_neigh_recvBuffer, h_migrating_sendBuffer, h_migrating_recvBuffer, h_mass_sendBuffer, h_mass_recvBuffer, flag_migration, buff_sizes, nparts_rank, d_add_mem);

    h_pack_particles(d, local_pos, local_vel, local_f, local_m, mi, h_neigh_sendBuffer, h_migrating_sendBuffer, h_mass_sendBuffer, boxdim, r_cut, r_skin, flag_migration);   



    // The exchange here is written explicitly inside the main code since is done 
    // while evaluating the pair forces for particles only inside the rank's local
    // domain.
    MPI_Status buff_status[12];
    MPI_Request buff_requests[12];
    int n_requests = 0;

    if(p.neighbour_east >= 0)
    {
      MPI_Isend(&h_neigh_sendBuffer[0], d.neigh_e_send, MPI_Field, p.neighbour_east, 0, p.comm, &buff_requests[n_requests]);
      n_requests += 1;
      MPI_Irecv(&h_neigh_recvBuffer[0], d.neigh_e_recv, MPI_Field, p.neighbour_east, 0, p.comm, &buff_requests[n_requests]);
      n_requests += 1;
      if(flag_migration){
          MPI_Isend(&h_migrating_sendBuffer[0], 3*d.migrating_e_send, MPI_Field, p.neighbour_east, 1, p.comm, &buff_requests[n_requests]);
          n_requests += 1;
          MPI_Irecv(&h_migrating_recvBuffer[0], 3*d.migrating_e_recv, MPI_Field, p.neighbour_east, 1, p.comm, &buff_requests[n_requests]);
          n_requests += 1;
          MPI_Isend(&h_mass_sendBuffer[0], d.migrating_e_send, MPI_DOUBLE, p.neighbour_east, 1, p.comm, &buff_requests[n_requests]);
          n_requests += 1;
          MPI_Irecv(&h_mass_recvBuffer[0], d.migrating_e_recv, MPI_DOUBLE, p.neighbour_east, 1, p.comm, &buff_requests[n_requests]);
          n_requests += 1;
        }
    }
    if(p.neighbour_west >= 0)
    {
      MPI_Isend(&h_neigh_sendBuffer[d.neigh_e_send], d.neigh_w_send, MPI_Field, p.neighbour_west, 0, p.comm, &buff_requests[n_requests]);
      n_requests += 1;
      MPI_Irecv(&h_neigh_recvBuffer[d.neigh_e_recv], d.neigh_w_recv, MPI_Field, p.neighbour_west, 0, p.comm, &buff_requests[n_requests]);
      n_requests += 1;
      if(flag_migration){
          MPI_Isend(&h_migrating_sendBuffer[3*d.migrating_e_send], 3*d.migrating_w_send, MPI_Field, p.neighbour_west, 1, p.comm, &buff_requests[n_requests]);
          n_requests += 1;
          MPI_Irecv(&h_migrating_recvBuffer[3*d.migrating_e_recv], 3*d.migrating_w_recv, MPI_Field, p.neighbour_west, 1, p.comm, &buff_requests[n_requests]);
          n_requests += 1;
          MPI_Isend(&h_mass_sendBuffer[d.migrating_e_send], d.migrating_w_send, MPI_DOUBLE, p.neighbour_west, 1, p.comm, &buff_requests[n_requests]);
          n_requests += 1;
          MPI_Irecv(&h_mass_recvBuffer[d.migrating_e_recv], d.migrating_w_recv, MPI_DOUBLE, p.neighbour_west, 1, p.comm, &buff_requests[n_requests]);
          n_requests += 1;
        }
    }


    // evaluation of local domain forces
    h_local_LJForces(local_f, local_pos, en, sigma, epsilon, r_cut, r_skin, boxdim, d.npart);

    // wait for all ranks to recieve all neighbouring and migrating particles
    MPI_Waitall(n_requests, buff_requests, MPI_STATUS_IGNORE);


    // evaluation of neighbouring forces
    h_neigh_LJForces(local_f, local_pos, h_neigh_recvBuffer, en, sigma, epsilon, r_cut, r_skin, boxdim, d.npart, d.neigh_e_recv+d.neigh_w_recv);  
    // final middle velocity step
    h_vl_step(local_vel, local_f, local_m, v2, timestep_length, d.npart);  




    // Fill into each rank domain the incoming particles by removing
    // the migrated ones
    if(flag_migration){
      for(int i = 0; i < d.migrating_e_send+d.migrating_w_send; i++){
        swap(local_pos[mi[i]], local_pos[d.npart - i - 1]);
        swap(local_vel[mi[i]], local_vel[d.npart - i - 1]);
        swap(local_f[mi[i]], local_f[d.npart - i - 1]);
        swap(local_m[mi[i]], local_m[d.npart - i - 1]);
      }
      int k = -(d.migrating_e_send+d.migrating_w_send);
      for(int i = 0; i < d.migrating_e_recv + d.migrating_w_recv; i++){
          local_pos[d.npart + k] = h_migrating_recvBuffer[i]; 
          local_vel[d.npart + k] = h_migrating_recvBuffer[i + d.migrating_e_recv+d.migrating_w_recv]; 
          local_f[d.npart + k] = h_migrating_recvBuffer[i + 2*(d.migrating_e_recv+d.migrating_w_recv)]; 
          local_m[d.npart + k] = h_mass_recvBuffer[i];
          k++;
      }
    }

  // update total number of local particles
  d.npart = nparts_rank[p.x];

  // stop timer before printing function
  nTimeStop = get_time ();      

  // printing section

    if(t % phys_out_freq == 0){
      global_en = 0;
      global_k = 0;
      if(d.npart > 0){
        MPI_Reduce(&en, &global_en, 1, MPI_DOUBLE, MPI_SUM, 0, p.comm);
        MPI_Reduce(&v2, &global_k, 1, MPI_DOUBLE, MPI_SUM, 0, p.comm);
      }
      else{
        v2 = 0;
        en = 0;
        MPI_Reduce(&en, &global_en, 1, MPI_DOUBLE, MPI_SUM, 0, p.comm);
        MPI_Reduce(&v2, &global_k, 1, MPI_DOUBLE, MPI_SUM, 0, p.comm);
      }
      if(p.x == 0)
         phys_out<< t << setw(15) << double((global_en + global_k)/N) << setw(15) << double(global_en/N) << setw(15) << double(global_k/N) <<endl;
    }

    if(t % vtk_out_freq == 0){
      string ss = "vtk_output/rank"+to_string(p.x)+vtk_out_name_base+to_string(vtk)+".vtk";
      ofstream vtk_file;
      vtk_file.open(ss.c_str());
      print_vtk(local_pos, local_vel, d.npart, local_m, vtk_file);
      vtk++;
    }

    nTimeLost += get_time () - nTimeStop;

  }   
// // ---------------------- End Loop


// // ----------------------------------------- End simulation
double global_time = 0;
double nTotalTime = get_time () - nTimeStart - nTimeLost;
MPI_Reduce(&nTotalTime, &global_time, 1, MPI_DOUBLE, MPI_SUM, 0, p.comm);

if(p.x == 0){
printf ("total time elapsed: %g milliseconds\n", (global_time) / 1000.0);    //milliseconds
}

  MPI_Finalize();
  return 0;

}
