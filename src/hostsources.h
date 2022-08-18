#ifndef HOSTSOURCES_H
#define HOSTSOURCES_H


#include <iostream>
#include <fstream>
#include <string.h>
#include <algorithm>
#include <mpi.h>
#include "mpiconsts.h"
#include "consts.h"


using namespace std;

// -------------- single-rank main functions
double h_boundary_ds(double x1, double x2, double r_skin, float dim);
void h_mv_step(Field* & Pos, Field* Vel, Field* Forces, double* m, double dt, Box boxdim, int N);
void h_vl_step(Field* & Vel, Field* Forces, double* m, double &v2, double dt, int N);
void h_local_LJForces(Field* & Forces, Field* Pos, double &en, double sigma, double epsilon, double r_cut, double r_skin, Box boxdim, int N);
void h_neigh_LJForces(Field* & Forces, Field* Pos, Field* Npos, double &en, double sigma, double epsilon, double r_cut, double r_skin, Box boxdim, int N, int N_neigh);
void h_LJForces(Field* & Forces, Field* Pos, Field* Npos, double &en, double sigma, double epsilon, double r_cut, double r_skin, Box boxdim, int N, int N_neigh);
void h_addLJForce(Field& f, Field pa, Field pb, double &en, double epsilon, double sigma, double r_cut, double r_skin, Box boxdim);
void h_add2LJForce(Field& fa, Field& fb, Field pa, Field pb, double &en, double epsilon, double sigma, double r_cut, double r_skin, Box boxdim);
void h_check_dims(Domain& d, Partition p, Field* &h_pos, Field* &h_vel, Field* &h_forces, double* &h_m, int* &mi, Field* &h_neigh_sendBuffer, Field* &h_neigh_recvBuffer, Field* &h_migrating_sendBuffer, Field* &h_migrating_recvBuffer, double* &h_mass_sendBuffer, double* &h_mass_recvBuffer, bool flag_migration, Sizes& buff_sizes, int* nparts_counter, float d_add_mem);
void h_pack_migrating(Domain d, Field* h_pos, Field* h_vel, Field* h_forces, double* h_m, Field* &h_migrating_sendBuffer, double* &h_migrating_mass_sendBuffer, Box boxdim, double r_cut, double r_skin);
void h_pack_neighbours(Domain d, Field* h_pos, Field* & h_neigh_sendBuffer, Sizes buff_sizes, double r_cut, double r_skin, Box boxdim);
void h_count_neighbours(Domain& d, Field* d_local_pos, Box boxdim, double r_cut, double r_skin);
void h_pack_particles(Domain d, Field* h_pos, Field* h_vel, Field* h_forces, double* h_m, int* &mi, Field* &h_neigh_sendBuffer, Field* &h_migrating_sendBuffer, double* &h_mass_sendBuffer, Box boxdim, double r_cut, double r_skin, bool flag_migration);

// Init functions
void loadparameters(string file_par, string* parameters);
void readparameters(string* parameters, string & input_file, double & time_end, \
  double & timestep_length,double & epsilon,double & sigma,int & part_out_freq, \
  string & part_out_name_base,int & vtk_out_freq, string & vtk_out_name_base,\
  int & cl_workgroup_1dsize, int & cl_workgroup_3dsize_x, \
  int & cl_workgroup_3dsize_y, int & cl_workgroup_3dsize_z,\
  float & x_min, float & x_max, float & y_min,\
  float & y_max, float & z_min, float & z_max,\
  int & x_n, int & y_n, int & z_n,\
  double & r_cut, double & r_skin);
void get_input_particles(Domain &d, string input_file, int& N, Field* &local_pos, Field* &local_vel, double* &local_m, float x_min, float y_min, float z_min);


//  Write functions
void print_phys(ofstream phys_out, double &en, double v2, double t);
void print_vtk(Field* Pos, Field* Vel, int N, double* m, ofstream & vtk_file);
void print_out(Field* Pos, Field* Vel, int N, double* m, ofstream & out_file);

void h_reallocateFieldCPU(Field* &ptr, int new_size);
void h_reallocateDoubleCPU(double* &ptr, int new_size);

// -------------- multi-rank main functions
void get_neighbours_info(Domain &d, Partition p, bool flag_migration);
void sendrecv_neighbours_info(Domain &d, Partition p, MPI_Request* request, int &n_requests, bool flag_migration);
void get_particles(Domain d, Partition p, Field* local_pos, Field* &neigh_sendBuff, Field* &neigh_recvBuff, Field* &migrating_sendBuff, Field* &migrating_recvBuff, double* &mass_sendBuff, double* &mass_recvBuff, MPI_Datatype MPI_Field, bool flag_migration);
void sendrecv_particles(Domain d, Partition p, Field* local_pos, Field* &neigh_sendBuff, Field* &neigh_recvBuff, Field* &migrating_sendBuff, Field* &migrating_recvBuff, double* &mass_sendBuff, double* &mass_recvBuff, MPI_Datatype MPI_Field, bool flag_migration, MPI_Status* status, MPI_Request* request, int& n_requests);


void get_migrating_particles(Domain d, Partition p, Field* &local_pos, Field* &local_vel, Field* &local_forces, Field* &migrating_sendBuff, Field* &migrating_recvBuff, MPI_Datatype MPI_Field);
void sendrecv_migrating_particles(Domain d, Partition p, Field* &migration_recvBuffer, Field* &migration_sendBuffer, MPI_Datatype MPI_Field, MPI_Request* request, int &n_requests);






#endif
