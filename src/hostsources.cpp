#include <cmath>
#include <iostream>
#include <fstream>
#include <string.h>
#include <algorithm>
#include <mpi.h>
#include "hostsources.h"
#define MPI_CHECK(call) \
    if((call) != MPI_SUCCESS) { \
        fprintf(stderr,"MPI error calling \""#call"\"\n");      \
        my_abort(-1); }

void my_abort(int err) {
      printf("Test FAILED\n");
      MPI_Abort(MPI_COMM_WORLD, err);
}

using namespace std;

Partition createPartition(int mpi_rank, int mpi_size)
{
    Partition p;

    // Determine size of the grid of MPI processes p.nx
    int dims[1] = {0};
    MPI_Dims_create(mpi_size, 1, dims);
    p.nx = dims[0];
    
    // Create cartesian communicator (p.comm), we do not allow the reordering of ranks here
    int periods[1] = {1};
    MPI_Comm comm_cart;
    MPI_Cart_create(MPI_COMM_WORLD, 1,  dims, periods, 0, &comm_cart);
    p.comm = comm_cart;
    
    // Determine the coordinate p.x
    int coordinate[1];
    MPI_Cart_coords(p.comm, mpi_rank, 1, coordinate);
    p.x = coordinate[0];
    
    // Determine neighbours
    // MPI_Cart_shift(p.comm, 0, 1, p.neighbour_south, p.neighbour_north);
    MPI_Cart_shift(p.comm, 1, 1, &p.neighbour_west, &p.neighbour_east);


    return p;
}


Domain createDomain(Partition p, Box boxdims)
{
    Domain d;
    
    // Size of the local domain
    if(boxdims.x / p.nx <= boxdims.r_min){
        printf("The input file chosen is meant to be run with less than %i processes.\n", round(boxdims.x /  boxdims.r_min));
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    d.nx = boxdims.x / p.nx;
    // Local domain left boundary
    d.startx = d.nx * p.x;
    // Local domain right boundary
    d.endx = d.startx + d.nx;
    
    return d;
}


// ------------------------------------- Host functions


// ----------- Linked List (disabled)

// void host_reset_cell_list(struct linked_list ** list,int N){

//   for(int i = 0; i <N; i++){
//     (*list)->cell_list[i] = -1;
//   }
// }


// void host_update_cell_list(struct linked_list ** list, Field* Pos, int N, Box boxdim){
  
//   unsigned m, n, k, cell_index;
//   for (int j = 0; j < N; j++){
//     m = floor(Pos[j].x/boxdim.x*boxdim.x_n);
//     n = floor(Pos[j].y/boxdim.y*boxdim.y_n);
//     k = floor(Pos[j].z/boxdim.z*boxdim.z_n);
//     cell_index = m + boxdim.x_n * (n + boxdim.z_n * k);
    
//     (*list)->particle_list[j] = j;
//     swap((*list)->particle_list[j], (*list)->cell_list[cell_index]);
//   }
// }

// struct linked_list * host_generate_list(Field* Pos, int x_n, int y_n, int z_n, Box &boxdim, int N, double r_cut, double r_skin)
//     // float cell_side_x = boxdim.x/x_n;
//     // float cell_side_y = boxdim.x/y_n;
//     // float cell_side_z = boxdim.x/z_n;

//     struct linked_list * list = (struct linked_list*)malloc(sizeof(struct linked_list));
//     boxdim.x_n = x_n;
//     boxdim.y_n = y_n;
//     boxdim.z_n = z_n;
//     boxdim.n_of_cells = x_n*y_n*z_n;

//     list->cell_list = (int*) malloc(boxdim.n_of_cells * sizeof(int));
//     list->particle_list = (int*) malloc(N * sizeof(int));
    
//     return list;
// }


// ----------- Leap Frog scheme

double h_boundary_ds(double x_1, double x_2, double r_skin, float dim){
  double x_ij = x_1 - x_2;
  return x_ij - dim*round(x_ij/dim); 
}

// move step including first half velocity update

void h_mv_step(Field* & Pos, Field* Vel, Field* Forces, double* m, double dt, Box boxdim, int N){
  double new_x, new_y, new_z;
  for(int i = 0; i < N; i++){
    //update position
    new_x = Pos[i].x + Vel[i].x*dt + Forces[i].x*pow(dt,2)/(m[i]*2); 
    new_y = Pos[i].y + Vel[i].y*dt + Forces[i].y*pow(dt,2)/(m[i]*2);
    new_z = Pos[i].z + Vel[i].z*dt + Forces[i].z*pow(dt,2)/(m[i]*2);
    Pos[i].x = fmod(new_x + 2*boxdim.x, boxdim.x);
    Pos[i].y = fmod(new_y + 2*boxdim.y, boxdim.y);
    Pos[i].z = fmod(new_z + 2*boxdim.z, boxdim.z);

    // first half of velocity update 
    Vel[i].x += Forces[i].x*dt/(m[i]*2);
    Vel[i].y += Forces[i].y*dt/(m[i]*2);
    Vel[i].z += Forces[i].z*dt/(m[i]*2);
  }
}

void h_vl_step(Field* & Vel, Field* Forces, double* m, double &v2, double dt, int N){
  for(int i=0;i<N;i++){
    Vel[i].x += Forces[i].x*dt/(m[i]*2);
    Vel[i].y += Forces[i].y*dt/(m[i]*2);
    Vel[i].z += Forces[i].z*dt/(m[i]*2);
    v2 += (pow(Vel[i].x, 2) + pow(Vel[i].y, 2) + pow(Vel[i].z, 2))/(2*m[i]);
  }
}


// ----------- Forces
void h_add2LJForce(Field& Force_1, Field& Force_2, Field Pos_1, Field Pos_2, double &en, double epsilon, double sigma, double r_cut, double r_skin, Box boxdim){
  double x = h_boundary_ds(Pos_1.x, Pos_2.x, r_skin, boxdim.x);
  double y = h_boundary_ds(Pos_1.y, Pos_2.y, r_skin, boxdim.y);
  double z = h_boundary_ds(Pos_1.z, Pos_2.z, r_skin, boxdim.z);
  double r = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
  if(r < r_cut){
    double q = (double) sigma/r;
    double qi = pow(q,6);
    double f = (double) 24*epsilon*pow(r,-2)*qi*(2*qi- 1);
    Force_1.x += f*x; Force_1.y += f*y; Force_1.z += f*z;
    Force_2.x += -f*x; Force_2.y += -f*y; Force_2.z += -f*z;
    en += 4*qi*(qi - 1); // update total energy
  }
}
void h_addLJForce(Field& Force_1, Field Pos_1, Field Pos_2, double &en, double epsilon, double sigma, double r_cut, double r_skin, Box boxdim){
  double x = h_boundary_ds(Pos_1.x, Pos_2.x, r_skin, boxdim.x);
  double y = h_boundary_ds(Pos_1.y, Pos_2.y, r_skin, boxdim.y);
  double z = h_boundary_ds(Pos_1.z, Pos_2.z, r_skin, boxdim.z);
  double r = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
  if(r < r_cut){
    double q = (double) sigma/r;
    double qi = pow(q, 6);
    double f = (double) 24*epsilon*pow(r,-2)*qi*(2*qi - 1);
    Force_1.x += f*x; Force_1.y += f*y; Force_1.z += f*z;
    en += 2*qi*(qi - 1); //single particle update of total energy
  }
}


void h_LJForces(Field* & Forces, Field* Pos, Field* Pos_neigh, double &en, double sigma, double epsilon, double r_cut, double r_skin, Box boxdim, int N, int N_neigh){
  for(int i=0;i<N;i++){
      Forces[i].x = 0; Forces[i].y = 0; Forces[i].z = 0;
      for(int j=0;j<i;j++){
        h_add2LJForce(Forces[i], Forces[j], Pos[i], Pos[j], en, epsilon, sigma, r_cut, r_skin, boxdim);
      }
      for(int k=0;k<N_neigh;k++){
        h_addLJForce(Forces[i], Pos[i], Pos_neigh[k], en, epsilon, sigma, r_cut, r_skin, boxdim);
      }
  }
}

void h_local_LJForces(Field* & Forces, Field* Pos, double &en, double sigma, double epsilon, double r_cut, double r_skin, Box boxdim, int N){
  for(int i=0;i<N;i++){
  Forces[i].x = 0; Forces[i].y = 0; Forces[i].z = 0;
    for(int j=0;j<i;j++){
      h_add2LJForce(Forces[i], Forces[j], Pos[i], Pos[j], en, epsilon, sigma, r_cut, r_skin, boxdim);
    }
  }
}

void h_neigh_LJForces(Field* & Forces, Field* Pos, Field* Pos_neigh, double &en, double sigma, double epsilon, double r_cut, double r_skin, Box boxdim, int N, int N_neigh){
  for(int i=0;i<N;i++){
    for(int j=0;j<N_neigh;j++){
      h_addLJForce(Forces[i], Pos[i], Pos_neigh[j], en, epsilon, sigma, r_cut, r_skin, boxdim);
    }
  }
}


// Linked list force evaluation (disabled)

// void host_linklist_LJForces(Field* & Forces, Field* Pos, linked_list list, double sigma, double epsilon, double r_cut, double r_skin, Box boxdim, int N){
//   unsigned cx, cy, cz;
//   int j;
//   for (int i = 0; i < N; i++){
//     cx = floor(Pos[i].x/boxdim.x*boxdim.x_n);
//     cy = floor(Pos[i].y/boxdim.y*boxdim.y_n);
//     cz = floor(Pos[i].z/boxdim.z*boxdim.z_n);
//     Forces[i].x = 0; Forces[i].y = 0; Forces[i].z = 0;
//     for(int off_x = -1; off_x <= 1; off_x++) {
//         for(int off_y = -1; off_y <= 1; off_y++) {
//             for(int off_z = -1; off_z <= 1; off_z++) {
//                 int ncx = (cx + off_x + 2 * boxdim.x_n) % boxdim.x_n;
//                 int ncy = (cy + off_y + 2 * boxdim.y_n) % boxdim.y_n;
//                 int ncz = (cz + off_z + 2 * boxdim.z_n) % boxdim.z_n;
//                 int nb_cell_id = ncx + boxdim.x_n * (ncy + boxdim.z_n * ncz);
//                 j = list.cell_list[nb_cell_id];
//                 while(j != -1){
//                   if(j != i)
//                     host_addLJForce(Forces[i], Pos[i], Pos[j], epsilon, sigma, r_cut, r_skin, boxdim);
//                   j = list.particle_list[j];
//                 }
//             }
//         }
//     }
//   }

// }







// ------- Memory management functions

void h_reallocateFieldCPU(Field* &ptr, int new_size){
  if(!ptr){
   ptr = (Field*)calloc(new_size, sizeof(Field));
  }
  Field* tmp = (Field*)realloc(ptr, new_size * sizeof(Field));
  if(tmp == NULL)
  {
      // Case 3, clean up then terminate.
      free(ptr);
      exit(0);
  }
  else if(tmp == ptr)
  {
      // Case 1: They point to the same place, so technically we can get away with
      // doing nothing.
      // Just to be safe, I'll assign NULL to tmp to avoid a dangling pointer.
      tmp = NULL;
  }
  else
  {
      // Case 2: Now tmp is a different chunk of memory.
      ptr = tmp;
      tmp = NULL;
  }
}
void h_reallocateDoubleCPU(double* &ptr, int new_size){
  if(!ptr){
    ptr = (double*)calloc(new_size, sizeof(double));
  }
  double* tmp = (double*)realloc(ptr, new_size * sizeof(double));
  if(tmp == NULL)
  {
      // Case 3, clean up then terminate.
      free(ptr);
      exit(0);
  }
  else if(tmp == ptr)
  {
      // Case 1: They point to the same place, so technically we can get away with
      // doing nothing.
      // Just to be safe, I'll assign NULL to tmp to avoid a dangling pointer.
      tmp = NULL;
  }
  else
  {
      // Case 2: Now tmp is a different chunk of memory.
      ptr = tmp;
      tmp = NULL;
  }
}




void h_check_dims(Domain &d, Partition p, Field* &h_pos, Field* &h_vel, Field* &h_forces, double* &h_m, int* &mi, Field* &h_neigh_sendBuffer, Field* &h_neigh_recvBuffer, Field* &h_migrating_sendBuffer, Field* &h_migrating_recvBuffer, double* &h_mass_sendBuffer, double* &h_mass_recvBuffer, bool flag_migration, Sizes& buff_sizes, int* nparts_counter, float d_add_mem){

  nparts_counter[p.x] = d.npart;
  if(d.neigh_e_send+d.neigh_w_send > buff_sizes.size_neigh_send){
    //reallocate memory for bigger buffer size
    int new_size = int(d_add_mem*(d.neigh_e_send+d.neigh_w_send));
    free(h_neigh_sendBuffer);
    h_neigh_sendBuffer = (Field*)malloc(new_size*sizeof(Field));
    //save new send buffer size
    buff_sizes.size_neigh_send = new_size;
  }
  if(d.neigh_e_recv+d.neigh_w_recv > buff_sizes.size_neigh_recv){
    //reallocate more memory for the incoming future neighbours
    int new_size = int(d_add_mem*(d.neigh_e_recv+d.neigh_w_recv));
    free(h_neigh_recvBuffer);
    h_neigh_recvBuffer = (Field*)malloc(new_size*sizeof(Field));
    //save new recv buffer size 
    buff_sizes.size_neigh_recv = new_size;
  }

  if(flag_migration){
    int r = d.migrating_w_recv+d.migrating_e_recv;
    int s = d.migrating_e_send+d.migrating_w_send;
    if(r > s)
    {
      if(d.npart + r - s > d.memory_size){
        int new_size = d.npart + int((r - s)*(d_add_mem));
        h_reallocateFieldCPU(h_pos, new_size);
        h_reallocateFieldCPU(h_vel, new_size);
        h_reallocateFieldCPU(h_forces, new_size);
        h_reallocateDoubleCPU(h_m, new_size);
        d.vacant = new_size - d.memory_size;
        nparts_counter[p.x] = d.npart + r - s;
        d.memory_size = new_size;
      }
      else{
        nparts_counter[p.x] = d.npart + r - s;
        d.vacant = d.memory_size - nparts_counter[p.x];
      }
    }
    if(s > r){
      nparts_counter[p.x] = d.npart - s + r;
      d.vacant = d.memory_size - nparts_counter[p.x];
    }
    if(s > buff_sizes.size_migrating_send)
    {
    //free past send_buffer pointers, no need to copy old data
      int new_size = int(d_add_mem*s);
      free(h_migrating_sendBuffer);
      free(h_mass_sendBuffer);
      free(mi);
      h_migrating_sendBuffer = (Field*)malloc(3*new_size*sizeof(Field));
      h_mass_sendBuffer = (double*)malloc(new_size*sizeof(double));
      mi = (int*) malloc (new_size*sizeof(int));
      //save new send_buffer size
      buff_sizes.size_migrating_send = new_size;
    }
    if(r > buff_sizes.size_migrating_recv)
    {
      //free past send_buffer pointers, no need to copy old data
      int new_size = int(d_add_mem*r);
      free(h_migrating_recvBuffer);
      h_migrating_recvBuffer = (Field*)malloc(3*new_size*sizeof(Field));
      h_mass_recvBuffer = (double*)malloc(new_size*sizeof(double));

      buff_sizes.size_migrating_recv = new_size;
    }
  }
  MPI_Allgather(&nparts_counter[p.x], 1, MPI_INT, nparts_counter, 1, MPI_INT, p.comm);
}


// ----------- Reading functions


// load parameters from input file
void loadparameters(string file_par, string* parameters){
  ifstream iFile;
  iFile.open(file_par.c_str(),ios::in);
  if (iFile.fail())
    exit(EXIT_FAILURE);

  string par_name;
  string value;
  int i = 0;
  while (iFile >> par_name >> value || i < 24)
    {
      parameters[i]=value;
      i++;
    }
  iFile.close();
}

// store input parameters to each variable
void readparameters(string* parameters, string & input_file, double & time_end, double & timestep_length,double & epsilon,double & sigma,int & part_out_freq, string & part_out_name_base,int & vtk_out_freq, string & vtk_out_name_base, int & cl_workgroup_1dsize, int & cl_workgroup_3dsize_x, int & cl_workgroup_3dsize_y, int & cl_workgroup_3dsize_z, float & x_min, float & x_max, float & y_min, float & y_max, float & z_min, float & z_max, int & x_n, int & y_n, int & z_n, double & r_cut, double & r_skin){
  
  input_file = parameters[0];
  time_end = stof(parameters[1]);
  timestep_length =stof(parameters[2]);

  epsilon =stof(parameters[3]);
  sigma=stof(parameters[4]);

  part_out_freq=stoi(parameters[5]);
  part_out_name_base=parameters[6];
  vtk_out_freq =stoi(parameters[7]);
  vtk_out_name_base= parameters[8];

  cl_workgroup_1dsize=stoi(parameters[9]);
  cl_workgroup_3dsize_x=stoi(parameters[10]);
  cl_workgroup_3dsize_y=stoi(parameters[11]);
  cl_workgroup_3dsize_z=stoi(parameters[12]);

  x_min=stof(parameters[13]);
  x_max=stof(parameters[14]);
  y_min=stof(parameters[15]);
  y_max=stof(parameters[16]);
  z_min=stof(parameters[17]);
  z_max=stof(parameters[18]);

  x_n=stoi(parameters[19]);
  y_n=stoi(parameters[20]);
  z_n=stoi(parameters[21]); 

  r_cut=stof(parameters[22]);
  r_skin=stof(parameters[23]);

}
 
// Read input file and assign particles to each rank
void get_input_particles(Domain& d, string input_file, int& N, Field* &local_pos, Field* &local_vel, double* &local_m, float x_min, float y_min, float z_min){

  input_file = "input/" + input_file;
  ifstream iFile;
  iFile.open(input_file.c_str(), ios::in);
  if (iFile.fail())
    exit(EXIT_FAILURE);
  
  iFile >> N;
  Field* all_pos = (Field*)calloc(N, sizeof(Field));
  Field* all_vel = (Field*)calloc(N, sizeof(Field));
  double* all_m = (double*)calloc(N, sizeof(double));

  int local_p = 0;
  int idx = 0;
  while (iFile >> all_m[idx] >> all_pos[idx].x >> all_pos[idx].y >> all_pos[idx].z >> \
    all_vel[idx].x >> all_vel[idx].y  >> all_vel[idx].z){
    // traslation of points into the box  B = [0, x_max - x_min] U [0, y_max - y_min] U ..
    all_pos[idx].x -= x_min;
    all_pos[idx].y -= y_min;
    all_pos[idx].z -= z_min;

    if(all_pos[idx].x < d.endx && all_pos[idx].x  >= d.startx){
      local_p++;
    }

    idx++;
  }

  d.npart = local_p;

  local_pos = (Field*)calloc(d.npart, sizeof(Field));
  local_vel = (Field*)calloc(d.npart, sizeof(Field));
  local_m = (double*)calloc(d.npart, sizeof(double));
  local_p = 0;
  for(int i = 0; i < N; i++){
    if(all_pos[i].x < d.endx && all_pos[i].x >= d.startx){
      local_pos[local_p] =  all_pos[i];
      local_vel[local_p] =  all_vel[i];
      local_m[local_p] =  all_m[i];
      local_p++;
    }
  }
    
  iFile.close();
  
}





// Functions for exchanging information

  void h_count_neighbours(Domain& d, Field* local_pos, Box boxdim, double r_cut, double r_skin){
    int neigh_east = 0;
    int neigh_west = 0;
    int migrating_e = 0;
    int migrating_w = 0;
    double d1, d2;
    for(int i = 0; i < d.npart; i++){
      d1 = h_boundary_ds(local_pos[i].x, d.endx, r_skin, boxdim.x);
      d2 = h_boundary_ds(local_pos[i].x, d.startx, r_skin, boxdim.x);
      if(abs(d1) <= r_cut)
        neigh_east++;
      if(d1 > 0)
        migrating_e++;
      if(abs(d2) <= r_cut)
        neigh_west++;
      if(d2 < 0)
        migrating_w++;
    }
    d.neigh_w_send = neigh_west;
    d.neigh_e_send = neigh_east;
    d.migrating_e_send = migrating_e;
    d.migrating_w_send = migrating_w;
  }



  void h_pack_particles(Domain d, Field* h_pos, Field* h_vel, Field* h_forces, double* h_m, int* &mi, Field* &h_neigh_sendBuffer, Field* &h_migrating_sendBuffer, double* &h_mass_sendBuffer, Box boxdim, double r_cut, double r_skin, bool flag_migration){
    int neigh_east = 0;
    int neigh_west = 0;
    int migrating_e = 0;
    int migrating_w = 0;
    double d1, d2;
    for(int i = 0; i < d.npart; i++){
      d1 = h_boundary_ds(h_pos[i].x, d.endx, r_skin, boxdim.x);
      d2 = h_boundary_ds(h_pos[i].x, d.startx, r_skin, boxdim.x);
      if(abs(d1) <= r_cut){
        h_neigh_sendBuffer[neigh_east]=h_pos[i];
        neigh_east++;
      }
      if(abs(d2) <= r_cut){
        h_neigh_sendBuffer[d.neigh_e_send + neigh_west]=h_pos[i];
        neigh_west++;
      }
      if(flag_migration){
        if(d1 > 0){
          mi[migrating_e]=i;
          h_migrating_sendBuffer[migrating_e]=h_pos[i];
          h_migrating_sendBuffer[migrating_e + d.migrating_e_send]=h_vel[i];
          h_migrating_sendBuffer[migrating_e + 2*d.migrating_e_send]=h_forces[i];
          h_mass_sendBuffer[migrating_e]=h_m[i];
          migrating_e++;
        }
        if(d2 < 0){
          mi[d.migrating_e_send + migrating_w]=i;
          h_migrating_sendBuffer[migrating_w + 3*d.migrating_e_send]=h_pos[i];
          h_migrating_sendBuffer[migrating_w + 3*d.migrating_e_send + d.migrating_w_send]=h_vel[i];
          h_migrating_sendBuffer[migrating_w + 3*d.migrating_e_send + 2*d.migrating_w_send]=h_forces[i];
          h_mass_sendBuffer[migrating_w + d.migrating_e_send]=h_m[i];
          migrating_w++;
        }
      }
    }
  }
  void get_neighbours_info(Domain &d, Partition p, bool flag_migration){

      MPI_Status status[8];
      MPI_Request request[8];
      int n_requests = 0;  
      sendrecv_neighbours_info(d, p,
        request, n_requests, flag_migration);

      MPI_CHECK(MPI_Waitall(n_requests, request, status));
       
  }
  void sendrecv_neighbours_info(Domain &d, Partition p,\
     MPI_Request* request, int &n_requests, bool flag_migration){
    
    if(p.neighbour_east >= 0){
      MPI_CHECK(MPI_Isend(&d.neigh_e_send, 1, MPI_INT, p.neighbour_east, 0, p.comm, &request[n_requests]));
      n_requests += 1;
      MPI_CHECK(MPI_Irecv(&d.neigh_e_recv, 1, MPI_INT, p.neighbour_east, 0, p.comm, &request[n_requests]));
      n_requests += 1;
      if(flag_migration){
        MPI_Isend(&d.migrating_e_send, 1, MPI_INT, p.neighbour_east, 0, p.comm, &request[n_requests]);
        n_requests += 1;
        MPI_Irecv(&d.migrating_e_recv, 1, MPI_INT, p.neighbour_east, 0, p.comm, &request[n_requests]);
        n_requests += 1;
      }
    }
    if(p.neighbour_west >= 0){
      MPI_CHECK(MPI_Isend(&d.neigh_w_send, 1, MPI_INT, p.neighbour_west, 0, p.comm, &request[n_requests]));
      n_requests += 1;
      MPI_CHECK(MPI_Irecv(&d.neigh_w_recv, 1, MPI_INT, p.neighbour_west, 0, p.comm, &request[n_requests]));
      n_requests += 1;
      if(flag_migration){
          MPI_Isend(&d.migrating_w_send, 1, MPI_INT, p.neighbour_west, 0, p.comm, &request[n_requests]);
          n_requests += 1;
          MPI_Irecv(&d.migrating_w_recv, 1, MPI_INT, p.neighbour_west, 0, p.comm, &request[n_requests]);
          n_requests += 1;
      }
    }

  }

  void h_pack_neighbours(Domain d, Field* h_pos,\
  Field* &h_neigh_sendBuffer, Sizes buff_sizes,
  double r_cut, double r_skin,\
  Box boxdim){

  //REMARKS:  
  //the selection includes also particles that already left the domain but are still inside the ghost 
  //strip given by r_min, so they will be included as neigh_pos for the correspoinding neighbouring
  //ranks. This allows the particle migration not to be computed every time step.
    int neigh_east = 0;
    int neigh_west = 0;
      for(int i = 0; i < d.npart; i++){
        if(abs(h_boundary_ds(h_pos[i].x, d.endx, r_skin, boxdim.x)) <= r_cut){
           h_neigh_sendBuffer[neigh_east]=h_pos[i];
            neigh_east++;
        }
        if(abs(h_boundary_ds(h_pos[i].x, d.startx, r_skin, boxdim.x)) <= r_cut){
            h_neigh_sendBuffer[d.neigh_e_send + neigh_west]=h_pos[i];
            neigh_west++;
        }
      }
  }

  void h_pack_migrating(Domain d,Field* h_pos, Field* h_vel, Field* h_forces, double* h_m, Field* &h_migrating_sendBuffer, double* &h_mass_sendBuffer, Box boxdim, double r_cut, double r_skin){
    int migrating_e = 0;
    int migrating_w = 0;
    for(int i = 0; i < d.npart; i++){
      if(h_boundary_ds(h_pos[i].x, d.endx, r_skin, boxdim.x) > 0)
      {
          h_migrating_sendBuffer[migrating_e]=h_pos[i];
          h_migrating_sendBuffer[migrating_e + d.migrating_e_send]=h_vel[i];
          h_migrating_sendBuffer[migrating_e + 2*d.migrating_e_send]=h_forces[i];
          h_mass_sendBuffer[migrating_e]=h_m[i];
          migrating_e++;
      }
      if(h_boundary_ds(h_pos[i].x, d.startx, r_skin, boxdim.x) < 0)
      {
        h_migrating_sendBuffer[migrating_w + 3*d.migrating_e_send]=h_pos[i];
        h_migrating_sendBuffer[migrating_w + 3*d.migrating_e_send + d.migrating_w_send]=h_vel[i];
        h_migrating_sendBuffer[migrating_w + 3*d.migrating_e_send + 2*d.migrating_w_send]=h_forces[i];
        h_mass_sendBuffer[migrating_w + d.migrating_e_send]=h_m[i];
        migrating_w++;
      }
    }
  }



  void get_particles(Domain d, Partition p, Field* local_pos, Field* &neigh_sendBuff,\
  Field* &neigh_recvBuff, Field* &migrating_sendBuff,\
  Field* &migrating_recvBuff, double* &mass_sendBuff,\
  double* &mass_recvBuff, MPI_Datatype MPI_Field, bool flag_migration){

    MPI_Status status[12];
    MPI_Request request[12];
    int n_requests = 0;  

    sendrecv_particles(d, p, local_pos,\
      neigh_sendBuff, neigh_recvBuff, migrating_sendBuff,\
      migrating_recvBuff, mass_sendBuff,\
      mass_recvBuff, MPI_Field, flag_migration, status, request, n_requests);

    MPI_Waitall(n_requests, request, status);
  }
    // Send the position of neighbouring particles.
  void sendrecv_particles(Domain d, Partition p, Field* local_pos, Field* &neigh_sendBuff,\
    Field* &neigh_recvBuff, Field* &migrating_sendBuff,\
    Field* &migrating_recvBuff, double* &mass_sendBuff,\
    double* &mass_recvBuff, MPI_Datatype MPI_Field, bool flag_migration, MPI_Status* status, MPI_Request* request, int& n_requests){

    if(p.neighbour_east >= 0) {
        MPI_Isend(&neigh_sendBuff[0], d.neigh_e_send, MPI_Field, p.neighbour_east, 0, p.comm, &request[n_requests]);
        n_requests += 1;
        MPI_Irecv(&neigh_recvBuff[0], d.neigh_e_recv, MPI_Field, p.neighbour_east, 0, p.comm, &request[n_requests]);
        n_requests += 1;
        if(flag_migration){
          MPI_Isend(&migrating_sendBuff[0], 3*d.migrating_e_send, MPI_Field, p.neighbour_east, 0, p.comm, &request[n_requests]);
          n_requests += 1;
          MPI_Irecv(&migrating_recvBuff[0], 3*d.migrating_e_recv, MPI_Field, p.neighbour_east, 0, p.comm, &request[n_requests]);
          n_requests += 1;
          MPI_Isend(&mass_sendBuff[0], d.migrating_e_send, MPI_DOUBLE, p.neighbour_east, 0, p.comm, &request[n_requests]);
          n_requests += 1;
          MPI_Irecv(&mass_recvBuff[0], d.migrating_e_recv, MPI_DOUBLE, p.neighbour_east, 0, p.comm, &request[n_requests]);
          n_requests += 1;
        }
    }
    
    
    if(p.neighbour_west >= 0){
        MPI_Isend(&neigh_sendBuff[d.neigh_e_send], d.neigh_w_send, MPI_Field,\
        p.neighbour_west, 0, p.comm, &request[n_requests]);
        n_requests += 1;
        MPI_Irecv(&neigh_recvBuff[d.neigh_e_recv], d.neigh_w_recv, MPI_Field,\
        p.neighbour_west, 0, p.comm, &request[n_requests]);
        n_requests += 1;
        if(flag_migration){
          MPI_Isend(&migrating_sendBuff[3*d.migrating_e_send], 3*d.migrating_w_send, MPI_Field, p.neighbour_west, 0, p.comm, &request[n_requests]);
          n_requests += 1;
          MPI_Irecv(&migrating_recvBuff[3*d.migrating_e_recv], 3*d.migrating_w_recv, MPI_Field, p.neighbour_west, 0, p.comm, &request[n_requests]);
          n_requests += 1;
          MPI_Isend(&mass_sendBuff[d.migrating_e_send], d.migrating_w_send, MPI_DOUBLE, p.neighbour_west, 0, p.comm, &request[n_requests]);
          n_requests += 1;
          MPI_Irecv(&mass_recvBuff[d.migrating_e_recv], d.migrating_w_recv, MPI_DOUBLE, p.neighbour_west, 0, p.comm, &request[n_requests]);
          n_requests += 1;
        }
    }

  }

// ----------- Print functions
// vtk unstructured grid
void print_vtk(Field* Pos, Field* Vel, int N, double* m, ofstream & vtk_file){
  vtk_file <<"# vtk DataFile Version 4.0 \nhesp visualization file \nASCII \nDATASET UNSTRUCTURED_GRID \nPOINTS "<<N<<" double \n";
  for (int i = 0; i < N; i++) {
    vtk_file <<fixed;
    vtk_file <<Pos[i].x<<" "<<Pos[i].y<<" "<<Pos[i].z<<endl;
  }
  vtk_file <<"CELLS "<<0<<" "<<0<<" \nCELL_TYPES "<<0<<"\nPOINT_DATA "<<N<<" \nSCALARS m double \nLOOKUP_TABLE default \n";
  for (int i = 0; i < N; i++) {
    vtk_file <<fixed;
    vtk_file <<m[i]<<endl;
  }
  vtk_file <<"VECTORS v double \n";
  for (int i = 0; i < N; i++) {
    vtk_file <<fixed;
    vtk_file <<Vel[i].x<<" "<<Vel[i].y<<" "<<Vel[i].z<<endl;
  }
};

void print_out(Field* Pos, Field* Vel, int N, double* m, ofstream & out_file){
  out_file <<N<<endl;
  for (int i=0; i<N; i++) out_file <<m[i]<<" "<<Pos[i].x<<" "<<Pos[i].y<< \
  " "<<Pos[i].z<<" "<<Vel[i].x<<" "<<Vel[i].y<<" "<<Vel[i].z<<endl;
};