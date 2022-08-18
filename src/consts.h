#ifndef CONSTS_H
#define CONSTS_H


struct Box
{
  float x; // boxdim in x direction
  float y; // boxdim in y direction
  float z;
  float volume;
  int n_of_cells;
  int x_n, y_n, z_n;
  double r_min;
};


struct Field
{
  double x;
  double y;
  double z;
};

typedef struct
{
    float nx;
    float startx;
    float endx;
    int npart;
    int vacant;
    int neigh_e_send, neigh_e_recv;
    int neigh_w_send, neigh_w_recv;
    int migrating_w_recv, migrating_w_send;
    int migrating_e_recv, migrating_e_send;
    int memory_size;
} Domain;

typedef struct
{
    int size_neigh_send;
    int size_neigh_recv;
    int size_migrating_send;
    int size_migrating_recv;
    int local_size;
} Sizes;



class Random {

	public:
		Random(int seed){
	    _seed = seed;
    }

		void SetA(int a) { _a = a;}
		void SetC(int c) { _c = c;}
		void SetM(int m) { _m = m;}
		
		double Rand(){
			double result = (_a*_seed + _c) % _m;
			_seed = result;
	
			return result/(_m-1);    //division returns a value between 0 e 1	
		}

	private:
	
		unsigned int _a;
		unsigned int _c;
		unsigned int _m;
		unsigned int _seed;

};



#endif
