#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <algorithm>
#include "consts.h"

using namespace std;



int main(int argc, char* argv[]){


  if(argc < 3){
		cout<<"Usage: ./initCubicGrid.o {npart density*} \n \n";
    abort();
  }

  int npart = atoi(argv[1]);
  double density = atof(argv[2]);



    // initilize output file
  ofstream Filein;
  Filein.open("input/generated_grid.in");
  if (Filein.fail())
    exit(EXIT_FAILURE);
  ofstream Filepar;
  Filepar.open("input/generated_grid.par");
  if (Filepar.fail())
    exit(EXIT_FAILURE);

    Filein << npart <<endl;
  

    // get box size
  float box_side = pow(npart/density,1./3);

    // find lowest perfect cube greater than or equal to the number of particles
  int nL = 2;
  while(pow(nL,3) < npart){
    nL = nL + 1;
  }

  // 3D index of the cubic lattice
  int index[3] = {0,0,0};

  // initialize Linear Congruential Generator for initial velocity distribution
  int seed = 3;
  Random LCG(seed);  
	LCG.SetA(1664525);
	LCG.SetC(1013904223);
	LCG.SetM(pow(2,31));

  Field* v = (Field*)malloc(npart*sizeof(Field));
  double sumv = 0;

    // writing to input file
    double vmax = 1000;
  for(int i = 0; i < npart; i++){
    v[i].x = vmax*(LCG.Rand() - 0.5); 
    v[i].y = vmax*(LCG.Rand() - 0.5); 
    v[i].z = vmax*(LCG.Rand() - 0.5); 
    sumv += v[i].x + v[i].y + v[i].z;
  }
  for(int i = 0; i < npart; i++){
    v[i].x = v[i].x/sumv; 
    v[i].y = v[i].y/sumv; 
    v[i].z = v[i].z/sumv; 
  }

  for (int i = 0; i < npart; i++)
  {
    Filein << 1 << " " << (index[0] + 0.5)*(box_side/nL) <<" "<< (index[1] + 0.5)*(box_side/nL) << \
    " " << (index[2] + 0.5)*(box_side/nL) << " " << v[i].x << " " << v[i].y << " " << v[i].z <<endl;


    // advance the index
    index[0] += 1;
    if(index[0] == nL){
      index[0] = 0;
      index[1] += 1;
      if(index[1] == nL){
        index[1] = 0;
        index[2] += 1;
      }
    }
  }



    // writing to parameter file
  Filepar << "part_input_file" << "  " << "generated_grid.in" << endl;
  Filepar << "time_end" << " " << 2 << endl;
  Filepar << "timestep_length" << " " << 0.001 << endl;
  Filepar << "epsilon" << " " << 1. << endl;
  Filepar << "sigma" << " " << 0.3 << endl;
  Filepar << "part_out_freq" << " " <<  10000000 << endl;
  Filepar << "part_out_name_base" << "  " << "generated_grid_" << endl;
  Filepar << "vtk_out_freq" << " " << 10 << endl;
  Filepar << "vtk_out_name_base" << " " << "generated_grid_" << endl;
  Filepar << "cl_workgroup_1dsize" << " " << 256 << endl;
  Filepar << "cl_workgroup_3dsize_x" << " " << 16 << endl;
  Filepar << "cl_workgroup_3dsize_y" << " " <<  16 << endl;
  Filepar << "cl_workgroup_3dsize_z" << " " << 1 << endl;
  Filepar << "x_min" << " " << 0 << endl;
  Filepar << "x_max" << " " <<  box_side << endl;
  Filepar << "y_min" << " " << 0 << endl;
  Filepar << "y_max" << " " <<  box_side << endl;
  Filepar << "z_min" << " " << 0 << endl;
  Filepar << "z_max" << " " <<  box_side << endl;
  Filepar << "x_n" << " " << 40 << endl;
  Filepar << "y_n" << " " << 40 << endl;
  Filepar << "z_n" << " " <<  3 << endl;
  Filepar << "r_cut" << " " <<  2.5 << endl;
  Filepar << "r_skin" << " " <<  0.3 << endl;


  return 0;

}