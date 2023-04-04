#include <stdio.h>
#include "common.h"


__global__ void  kernel_update_velocity(float* vx, float* vz, float* txx, float* tzz, float* txz,\
	float* C, float* pml_xtxx, float* pml_xtxz, float* pml_ztxz, float* pml_ztzz, float* pml_para,\
  int*  rec_loc, int* src_loc, float wave,\
	int NZ, int NX, int NK, int NR, int nabs, int fnabs,\
	float dx, float dz, float dt, int it){


  //printf("device vx[100] = %f\n", vx[100]);
  //vx[100] = 13;
  //printf("device vx[100] = %f\n", vx[100]);

	//local parameters
	int ix, iz, pos, posc, ipml, zcom, xcom, pmlpos;
	float xtxx, xtxz, ztxz, ztzz, Ttmp;
	float dtdx, dtdz;

	dtdx = dt/dx;
  dtdz = dt/dz;
	//add source
  // vz[ (src_loc[0] + nabs) * NZ + src_loc[1] + fnabs] += wave;
	// vx[ (src_loc[0] + nabs) * NZ + src_loc[1] + fnabs] += wave;
  
  ix = blockIdx.x * blockDim.x + threadIdx.x;
  iz = blockIdx.y * blockDim.y + threadIdx.y;
  //printf("%d\n", threadIdx.x);
  //printf("it %d\n", it);
  //printf("ix %d\n", ix);
  //printf("NX %d NZ %d in kernel_velocity\n", NX, NZ);
  if(ix > 0 &&  ix < NX - 1 && iz > 0 && iz < NZ - 1){
  //printf("ix %d\n", ix);
	// calculate derivative;
	//for(ix = 1; ix < NX - 1; ix++) 
	//    for(iz = 1; iz < NZ - 1; iz++) {
	    	pos = ix * NZ + iz;
	    	posc = ix * NZ * NK + iz * NK;

			if ( ix == 1 || ix == NX - 2 ) {
			xtxx = txx[pos] - txx[pos - NZ];
			xtxz = txz[pos + NZ] - txz[pos];
			}
			else {
			xtxx = coe1 * (txx[pos]  - txx[pos - NZ]) + coe2 * (txx[pos + NZ] - txx[pos - 2*NZ]);
			xtxz = coe1 * (txz[pos + NZ] - txz[pos])   + coe2 * (txz[pos + 2*NZ] - txz[pos - NZ]);
			}

			if( iz == 1 || iz == NZ - 2 )  {
			ztxz = txz[pos] - txz[pos - 1];
			ztzz = tzz[pos + 1] - tzz[pos];
			}
			else {
			ztxz = coe1 * (txz[pos]  - txz[pos - 1]) + coe2 * (txz[pos + 1] - txz[pos - 2]);
			ztzz=  coe1 * (tzz[pos + 1] - tzz[pos])  + coe2 * (tzz[pos + 2] - tzz[pos - 1]);
			}
			// dtdx = dt/dx
			vx[pos] += (xtxx * dtdx + ztxz * dtdz) * 2.0 / (C[posc + 0] + C[posc - NZ*NK + 0]);
			vz[pos] += (xtxz * dtdx + ztzz * dtdz) * 2.0 / (C[posc + 0] + C[posc + NK + 0]);
      
      //__syncthreads();

			zcom =  ((ix >= nabs) && (ix < NX - nabs) && (iz >= NZ - nabs)) ? 1:0;
			xcom =  (ix <= nabs) ? 0: ( (ix < NX - nabs) ? (ix - nabs) : (NX - 2*nabs));
			pmlpos = pos - ( NZ - fnabs - nabs) * (zcom + xcom);

		// PML at X direction
			if(ix < nabs) {
			 ipml = nabs - 1 - ix;
			 // vx
			 Ttmp = (2.0 * pml_xtxx[pmlpos] + pml_para[ipml*12+4] * xtxx) / (2.0 + pml_para[ipml*12+0]);
			 vx[pos] += ((pml_para[ipml*12+2] - 1) * xtxx / dx - pml_para[ipml*12+2] * Ttmp) * dt * 2.0/(C[posc + 0] + C[posc - NZ*NK + 0]);
			 pml_xtxx[pmlpos] = 2.0 * Ttmp - pml_xtxx[pmlpos];
			 // vz
			 Ttmp = (2.0 * pml_xtxz[pmlpos] + pml_para[ipml*12+5] * xtxz) / (2.0 + pml_para[ipml*12+1]);
			 //vz[pos] += ((pml_beta_half[ipml] - 1) * xtxz / dx - pml_beta_half[ipml] * Ttmp) * dt * 2.0/(C[posc + 0] + C[posc + NK + 0]);
			 vz[pos] += ((pml_para[ipml*12+3] - 1) * xtxz / dx - pml_para[ipml*12+3] * Ttmp) * dt * 2.0/(C[posc + 0] + C[posc + NK + 0]);
			 pml_xtxz[pmlpos] = 2.0 * Ttmp - pml_xtxz[pmlpos];
			}
			if(ix >= NX - nabs) {
			 // vz
			 ipml = ix - (NX - nabs);
			 Ttmp = (2.0 * pml_xtxz[pmlpos] + pml_para[ipml*12+5] * xtxz) / (2.0 + pml_para[ipml*12+1]);
			 //vz[pos] += ((pml_beta_half[ipml] - 1) * xtxz / dx - pml_beta_half[ipml] * Ttmp) * dt * 2.0/(C[posc + 0] + C[posc + NK + 0]);
			 vz[pos] += ((pml_para[ipml*12+3] - 1) * xtxz / dx - pml_para[ipml*12+3] * Ttmp) * dt * 2.0/(C[posc + 0] + C[posc + NK + 0]);
			 pml_xtxz[pmlpos] = 2.0 * Ttmp - pml_xtxz[pmlpos];
			}
			if(ix > NX - nabs) {
			 // vx
			 ipml = ix - (NX - nabs) - 1;
			 Ttmp = (2.0 * pml_xtxx[pmlpos] + pml_para[ipml*12+4] * xtxx) / (2.0 + pml_para[ipml*12+0]);
			 vx[pos] += ((pml_para[ipml*12+2] - 1) * xtxx / dx - pml_para[ipml*12+2] * Ttmp) * dt * 2.0/(C[posc + 0] + C[posc - NZ*NK + 0]);
			 pml_xtxx[pmlpos] = 2.0 * Ttmp - pml_xtxx[pmlpos];
			}

		// PML at Z direction
			if(iz < fnabs) {
			 ipml = fnabs - 1 - iz;
			 // vx
			 Ttmp = (2.0 * pml_ztxz[pmlpos] + pml_para[ipml*12+10] * ztxz) / (2.0 + pml_para[ipml*12+6]);
			 vx[pos] += ((pml_para[ipml*12+8] - 1) * ztxz / dz - pml_para[ipml*12+8] * Ttmp) * dt * 2.0/(C[posc + 0] + C[posc - NZ*NK + 0]);
			 pml_ztxz[pmlpos] = 2.0 * Ttmp - pml_ztxz[pmlpos];
			 // vz
			 Ttmp = (2.0 * pml_ztzz[pmlpos] + pml_para[ipml*12+11] * ztzz) / (2.0 + pml_para[ipml*12+7]);
			 //vz[pos] += ((pml_beta_half[ipml] - 1) * ztzz / dx - pml_beta_half[ipml] * Ttmp) * dt * 2.0/(C[posc + 0] + C[posc - NZ*NK + 0]);
			 vz[pos] += ((pml_para[ipml*12+9] - 1) * ztzz / dz - pml_para[ipml*12+9] * Ttmp) * dt * 2.0/(C[posc + 0] + C[posc - NZ*NK + 0]);
			 pml_ztzz[pmlpos] = 2.0 * Ttmp - pml_ztzz[pmlpos];
			}

			if(iz >= NZ - nabs) {
			 // vz
			 ipml = iz - (NZ - nabs);
			 Ttmp = (2.0 * pml_ztzz[pmlpos] + pml_para[ipml*12+11] * ztzz) / (2.0 + pml_para[ipml*12+7]);
			 //vz[pos] += ((pml_beta_half[ipml] - 1) * ztzz / dx - pml_beta_half[ipml] * Ttmp) * dt * 2.0/(C[posc + 0] + C[posc - NZ*NK + 0]);
			 vz[pos] += ((pml_para[ipml*12+9] - 1) * ztzz / dz - pml_para[ipml*12+9] * Ttmp) * dt * 2.0/(C[posc + 0] + C[posc - NZ*NK + 0]);
			 pml_ztzz[pmlpos] = 2.0 * Ttmp - pml_ztzz[pmlpos];
			}
			if(iz > NZ - nabs) {
			 // vx
			 ipml = iz - (NZ - nabs) - 1;
			 Ttmp = (2.0 * pml_ztxz[pmlpos] + pml_para[ipml*12+10] * ztxz) / (2.0 + pml_para[ipml*12+6]);
			 vx[pos] += ((pml_para[ipml*12+8] - 1) * ztxz / dz - pml_para[ipml*12+8] * Ttmp) * dt * 2.0/(C[posc + 0] + C[posc - NZ*NK + 0]);
			 pml_ztxz[pmlpos] = 2.0 * Ttmp - pml_ztxz[pmlpos];

			}

      // __syncthreads();
		}

	return;

}

/*
float* get_deviceMem_float(float* ori, int msize){

  float* dev_ori;
  cudaError_t err = cudaSuccess;
  err = cudaMalloc( (void**)&dev_ori, msize * sizeof(float) );
  if (err != cudaSuccess) {
    printf("stderr, %s\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }

  err = cudaMemcpy(dev_ori, ori, msize * sizeof(float), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) {
    printf("stderr, %s\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }

  return dev_ori;
}

int* get_deviceMem_int(int* ori, int msize){

  int* dev_ori;
  cudaError_t err = cudaSuccess;
  err = cudaMalloc( (void**)&dev_ori, msize * sizeof(int) );
  if (err != cudaSuccess) {
    printf("stderr, %s\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }

  err = cudaMemcpy(dev_ori, ori, msize * sizeof(int), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) {
    printf("stderr, %s\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }

  return dev_ori;
}

void copyback_float(float* ori, float* dev_ori, int msize){

  cudaError_t err = cudaSuccess;
  err = cudaMemcpy(ori, dev_ori, msize*sizeof(float), cudaMemcpyDeviceToHost);
  if (err != cudaSuccess) {
    printf("stderr, %s\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }

}
	
//void copyback_int(ori, dev_ori, msize){
//
//  err = cudaMemcpy(ori, dev_ori, msize*sizeof(int), cudaMemcpyDeviceToHost);
//  if (err != cudaSuccess) {
//    printf("stderr, %s\n", cudaGetErrorString(err));
//    exit(EXIT_FAILURE);
//  }
//}

	
void gpu_update_velocity(float* vx, float* vz, float* txx, float* tzz, float* txz,\
	float* C, float* pml_vxx, float* pml_vzz, float* pml_vxz, float* pml_vzx,\
	float* pml_xtxx, float* pml_xtxz, float* pml_ztxz, float* pml_ztzz, float* pml_para,\
	int NZ, int NX, int NK, int NR, int nabs, int fnabs,\
  int* rec_loc, float* record_vz, float* record_vx, float* div, float* curl,\
  float wave, int* src_loc,\
  int it, float dx, float dz, float dt, int pml_size){

  
  //printf("in host\n");
  cudaSetDevice(6);
  //cudaError_t err = cudaSuccess;
  int ir, pos;
  // local parameters
  float *dev_vx = get_deviceMem_float(vx, NX*NZ);
  float *dev_vz = get_deviceMem_float(vz, NX*NZ);
  float *dev_txx = get_deviceMem_float(txx, NX*NZ);
  float *dev_tzz = get_deviceMem_float(tzz, NX*NZ);
  float *dev_txz = get_deviceMem_float(txz, NX*NZ);
  float *dev_C = get_deviceMem_float(C, NX*NZ*NK);
  float *dev_pml_vxx = get_deviceMem_float(pml_vxx, pml_size);
  float *dev_pml_vzz = get_deviceMem_float(pml_vzz, pml_size);
  float *dev_pml_vxz = get_deviceMem_float(pml_vxz, pml_size);
  float *dev_pml_vzx = get_deviceMem_float(pml_vzx, pml_size);
  float *dev_pml_xtxx = get_deviceMem_float(pml_xtxx, pml_size);
  float *dev_pml_xtxz = get_deviceMem_float(pml_xtxz, pml_size);
  float *dev_pml_ztxz = get_deviceMem_float(pml_ztxz, pml_size);
  float *dev_pml_ztzz = get_deviceMem_float(pml_ztzz, pml_size);
  float *dev_pml_para = get_deviceMem_float(pml_para, pml_size);
  int *dev_rec_loc = get_deviceMem_int(rec_loc, NR);
  int *dev_src_loc = get_deviceMem_int(src_loc, 2);
  

  int block_size = 20;
  dim3 Threads(block_size, block_size, 1);
  dim3 Grids((NX + block_size - 1) / block_size, (NZ + block_size - 1) / block_size, 1);
  // cal
  //printf("in host dev_vx[100] = %f\n", dev_vx[100]);
 // kernel_test <<< 1,1 >>> (dev_vx);
  //kernel_update_velocity <<< Grids, Threads >>> ( dev_vx,  dev_vz,  dev_txx,  dev_tzz,  dev_txz,
  kernel_update_velocity <<< 1,  >>> ( dev_vx,  dev_vz,  dev_txx,  dev_tzz,  dev_txz,\
    dev_C, dev_pml_vxx,  dev_pml_vzz,  dev_pml_vxz,  dev_pml_vzx,\
	  dev_pml_xtxx,  dev_pml_xtxz,  dev_pml_ztxz,  dev_pml_ztzz,  dev_pml_para,\
    dev_rec_loc,  dev_src_loc,  wave,\
	  NZ,  NX,  NK,  NR,  nabs,  fnabs,\
	  dx,  dz,  dt);
  //printf("in host dev_vx[100] = %f\n", dev_vx[100]);

  // copy device parameters to host
  copyback_float(vx, dev_vx, NX*NZ);
  copyback_float(vz, dev_vz, NX*NZ);
  copyback_float(txx, dev_txx, NX*NZ);
  copyback_float(txz, dev_txz, NX*NZ);
  copyback_float(tzz, dev_tzz, NX*NZ);
  //copyback_float(C, dev_vx, NX*NZ*NK);
  copyback_float(pml_vxx, dev_pml_vxx, pml_size);
  copyback_float(pml_vzz, dev_pml_vzz, pml_size);
  copyback_float(pml_vxz, dev_pml_vxz, pml_size);
  copyback_float(pml_vzx, dev_pml_vzx, pml_size);
  copyback_float(pml_xtxx, dev_pml_xtxx, pml_size);
  copyback_float(pml_xtxz, dev_pml_xtxz, pml_size);
  copyback_float(pml_ztxz, dev_pml_ztzz, pml_size);
  copyback_float(pml_ztzz, dev_pml_ztzz, pml_size);
//  copyback_float(pml_para, dev_pml_para, pml_size);

//  printf("in host vx[100] = %f\n", vx[100]);

	// calculate the output data
	for(ir = 0; ir < NR; ir++){
		pos = (rec_loc[ir * 2 + 0] + nabs) * NZ + rec_loc[ir * 2 + 1] + fnabs;
//		record_vx[it * NR + ir] = vx[pos];
//		record_vz[it * NR + ir] = vz[pos];
		record_vx[it * NR + ir] = 0.5 * ( vx[pos] + vx[pos + NZ] );
		record_vz[it * NR + ir] = 0.5 * ( vz[pos] + vz[pos - 1] );
		div[it * NR + ir] =  (coe1 * (vx[pos + NZ] - vx[pos])
			 			    + coe2 * (vx[pos + 2*NZ] - vx[pos - NZ])) / dx
			 				+(coe1 * (vz[pos] - vz[pos - 1])
			 				+ coe2 * (vz[pos + 1] - vz[pos - 2])) / dz;
		curl[it * NR + ir] = 0.25 * (
			( (coe1*(vx[pos] - vx[pos-1])+coe2*(vx[pos+1] - vx[pos-2]))/dz) 
			- ( (coe1*(vz[pos-1] - vz[pos-NZ-1])+coe2*(vz[pos+NZ-1] - vz[pos-2*NZ-1]))/dx)
			+( (coe1*(vx[pos+1] - vx[pos])+coe2*(vx[pos+2] - vx[pos-1]))/dz) 
			- ( (coe1*(vz[pos] - vz[pos-NZ])+coe2*(vz[pos+NZ] - vz[pos-2*NZ]))/dx)
			+( (coe1*(vx[pos+NZ] - vx[pos+NZ-1])+coe2*(vx[pos+NZ+1] - vx[pos+NZ-2]))/dz) 
			- ( (coe1*(vz[pos+NZ-1] - vz[pos-1])+coe2*(vz[pos+2*NZ-1] - vz[pos-NZ-1]))/dx)
			+( (coe1*(vx[pos+NZ+1] - vx[pos+NZ])+coe2*(vx[pos+NZ+2] - vx[pos+NZ-1]))/dz) 
			- ( (coe1*(vz[pos+NZ] - vz[pos])+coe2*(vz[pos+2*NZ] - vz[pos-NZ]))/dx)
			);
  }

//  free(dev_vx);
  cudaFree(dev_vx);
  cudaFree(dev_vz);
  cudaFree(dev_txx);
  cudaFree(dev_txz);
  cudaFree(dev_tzz);
  cudaFree(dev_C);
  cudaFree(dev_pml_vxx);
  cudaFree(dev_pml_vzz);
  cudaFree(dev_pml_vxz);
  cudaFree(dev_pml_vzx);
  cudaFree(dev_pml_xtxx);
  cudaFree(dev_pml_xtxz);
  cudaFree(dev_pml_ztxz);
  cudaFree(dev_pml_ztzz);
  cudaFree(dev_pml_para);
  cudaFree(dev_rec_loc);
  cudaFree(dev_src_loc);

  }
*/
	
