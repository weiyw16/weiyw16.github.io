#define coe1      1.125               // 9/8 = 1.125
#define coe2      -0.041666666666666  // 1/24 = 0.041666666666666
#define pi        3.141592653589793
#define erro      0.0000000001


// functions


float* init_matrix(int sizenum);
int* init_matrix_int(int sizenum);
float Ricker(float time, float freq0, float amp0 );
void get_C(FILE* mdfid, int NK, int NZ, int NX, int nabs, int fnabs, int readin_model,\
    float* C, float& Vmax, float& Vmin );
void set_pml(int nabs, float dx, float dz, float dt, float Vmax, float freq0, float* pml_para);
void readin_coor(FILE* srfid, int NR, int* rec_loc, float dx, float dz);
void mymemorycp(float* buf_array, float* tar, int start, int len);



void update_velocity(float* vx, float* vz, float* txx, float* tzz, float* txz,\
	float* C, float* pml_vxx, float* pml_vzz, float* pml_vxz, float* pml_vzx,\
	float* pml_xtxx, float* pml_xtxz, float* pml_ztxz, float* pml_ztzz, float* pml_para,\
	int NZ, int NX, int NK, int NR, int nabs, int fnabs,\
  int* rec_loc, float* record_vz, float* record_vx, float* div, float* curl,\
  float wave, int* src_loc,\
	int it, float dx, float dz, float dt);
void update_stress(float* vx, float* vz, float* txx, float* tzz, float* txz,\
	float* C, float* pml_vxx, float* pml_vzz, float* pml_vxz, float* pml_vzx, float* pml_para,\
	int NZ, int NX, int NK, int nabs, int fnabs,\
	float wave, int* src_loc, float dx, float dz, float dt);
void gpu_update_velocity(float* vx, float* vz, float* txx, float* tzz, float* txz,\
	float* C, float* pml_vxx, float* pml_vzz, float* pml_vxz, float* pml_vzx,\
	float* pml_xtxx, float* pml_xtxz, float* pml_ztxz, float* pml_ztzz, float* pml_para,\
	int NZ, int NX, int NK, int NR, int nabs, int fnabs,\
  int* rec_loc, float* record_vz, float* record_vx, float* div, float* curl,\
  float wave, int* src_loc,\
	int it, float dx, float dz, float dt, int pml_size);
/*
__global__ void  kernel_update_velocity(float* vx, float* vz, float* txx, float* tzz, float* txz,\
  float* C, float* pml_vxx, float* pml_vzz, float* pml_vxz, float* pml_vzx,\
  float* pml_xtxx, float* pml_xtxz, float* pml_ztxz, float* pml_ztzz, float* pml_para,\
  int*  rec_loc, int* src_loc, float wave,\
  int NZ, int NX, int NK, int NR, int nabs, int fnabs,\
  float dx, float dz, float dt);

__global__ void kernel_update_stress (float* vx, float* vz, float* txx, float* tzz, float* txz,\
	float* C, float* pml_vxx, float* pml_vzz, float* pml_vxz, float* pml_vzx, float* pml_para,\
	int NZ, int NX, int NK, int nabs, int fnabs,\
	float wave, int* src_loc, float dx, float dz, float dt);
  */
