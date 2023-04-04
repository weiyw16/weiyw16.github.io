__global__ void  kernel_update_velocity(float* vx, float* vz, float* txx, float* tzz, float* txz,\
  float* C, float* pml_xtxx, float* pml_xtxz, float* pml_ztxz, float* pml_ztzz, float* pml_para,\
  int*  rec_loc, int* src_loc, float wave,\
  int NZ, int NX, int NK, int NR, int nabs, int fnabs,\
  float dx, float dz, float dt, int it);

__global__ void kernel_update_stress (float* vx, float* vz, float* txx, float* tzz, float* txz,\
	float* C, float* pml_vxx, float* pml_vzz, float* pml_vxz, float* pml_vzx, float* pml_para,\
	int NZ, int NX, int NK, int nabs, int fnabs,\
	float wave, int* src_loc, float dx, float dz, float dt);

//__global__ void kernel_update (float* vx, float* vz, float* txx, float* tzz, float* txz,\
//	float* C, float* pml_vxx, float* pml_vzz, float* pml_vxz, float* pml_vzx, \
//  float* pml_xtxx, float* pml_xtxz, float* pml_ztxz, float* pml_ztzz, float* pml_para,\
//  int*  rec_loc, int* src_loc, float wave,\
//	int NZ, int NX, int NK, int NR, int nabs, int fnabs,\
//	float dx, float dz, float dt);
