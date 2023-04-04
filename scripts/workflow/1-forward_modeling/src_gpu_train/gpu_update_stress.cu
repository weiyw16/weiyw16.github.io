#include <stdio.h>
#include "common.h"

/*
float get_cc(float av, float bv, float cv, float dv){

  float cc = 0;
  float erro = 0.0000000001;
  if ( av <= erro || bv <= erro || cv <= erro || dv <= erro )
    cc = 0.;
  else {  //c55=4.0/(1.0/av+1.0/bv+1.0/c+1.0/dv)
    cc = (4.0*av*bv*cv*dv ) / (bv*cv*dv + av*cv*dv + av*bv*dv + av*bv*cv);
  }
  return cc;
}
*/
__global__ void kernel_update_stress (float* vx, float* vz, float* txx, float* tzz, float* txz,\
	float* C, float* pml_vxx, float* pml_vzz, float* pml_vxz, float* pml_vzx, float* pml_para,\
	int NZ, int NX, int NK, int nabs, int fnabs,\
	float wave, int* src_loc, float dx, float dz, float dt){

	// printf("load in stress update\n");
	//local parameters
	int ix, iz, pos, posc, ipml, zcom, xcom, pmlpos;
	float Ctmp, Vtmp;
	float Dvxdz, Dvxdx, Dvzdz, Dvzdx;
  float dtdx, dtdz;

  dtdx = dt/dx;
  dtdz = dt/dz;

  //__syncthreads();

  ix = blockIdx.x * blockDim.x + threadIdx.x;
  iz = blockIdx.y * blockDim.y + threadIdx.y;

	//add source
  if (ix == ( src_loc[0] + nabs ) && iz == ( src_loc[1] + fnabs ) ){
	  txx[ (src_loc[0] + nabs) * NZ + src_loc[1] + fnabs] += wave;
	  tzz[ (src_loc[0] + nabs) * NZ + src_loc[1] + fnabs] += wave;
  }
  //printf("%d\n", blockIdx.x);
	// calculate stress
  if(ix > 0 &&  ix < NX - 1 && iz > 0 && iz < NZ - 1){
//	for (ix = 1; ix < NX - 1; ix++)
//		for (iz = 1; iz < NZ - 1; iz++){
			pos = ix * NZ + iz;
			posc = ix * NZ * NK + iz * NK;


			if ( ix == 1 || ix == NX - 2 ) {
			   Dvxdx = vx[pos+NZ] - vx[pos];//m!=mm-1
			   Dvzdx = vz[pos] - vz[pos-NZ];//m!=0
			}
			else {
			   Dvxdx = coe1*(vx[pos+NZ]-vx[pos]) + coe2*(vx[pos+2*NZ]-vx[pos-NZ]);  /* for Txx and Tzz, at (m+1/2,k) */
			   Dvzdx = coe1*(vz[pos]-vz[pos-NZ]) + coe2*(vz[pos+NZ]-vz[pos-2*NZ]);  /* for Txz, at (m,k+1/2) */
			}

			if( iz == 1 || iz == NZ - 2 ) {
			   Dvzdz = vz[pos] - vz[pos-1];//k!=0
			   Dvxdz = vx[pos+1] - vx[pos];//k!=kk-1
			}
			else {
			   Dvzdz = coe1*(vz[pos]-vz[pos-1]) + coe2*(vz[pos+1]-vz[pos-2]);   /* for Txx and Tzz, at (m+1/2,k) */
			   Dvxdz = coe1*(vx[pos+1]-vx[pos]) + coe2*(vx[pos+2]-vx[pos-1]);   /* for Txz, at (m,k+1/2) */
			}



//      Ctmp = get_cc( C[posc + 3], C[posc + 3 - NZ*NK], C[posc + NK + 3], C[posc + NK + 3 - NZ*NK]);
      if( C[posc + 3] <= erro || C[posc + 3 - NZ*NK] <= erro || C[posc + NK + 3] <= erro || C[posc + NK + 3 - NZ*NK] <= erro )
        Ctmp = 0;
      else {
        //Ctmp = 4.0 / ( 1.0 / C[posc + 3] + 1.0 / C[posc + 3 - NZ*NK] + 1.0 / C[posc + NK + 3] + 1.0 / C[posc + NK + 3 - NZ*NK])
        Ctmp = ( 4.0 * C[posc + 3] * C[posc + 3 - NZ*NK] * C[posc + NK + 3] * C[posc + NK + 3 - NZ*NK] ) \
               / ( C[posc + 3 - NZ*NK] * C[posc + NK + 3] * C[posc + NK + 3 - NZ*NK] \
                   + C[posc + 3] * C[posc + NK + 3] * C[posc + NK + 3 - NZ*NK] \
                   + C[posc + 3] * C[posc + 3 - NZ*NK] * C[posc + NK + 3 - NZ*NK] \
                   + C[posc + 3] * C[posc + 3 - NZ*NK] * C[posc + NK + 3] );
      }

			txx[pos] += C[posc + 1] * dtdx * Dvxdx + C[posc + 2] * dtdz * Dvzdz;
		  tzz[pos] += C[posc + 2] * dtdx * Dvxdx + C[posc + 1] * dtdz * Dvzdz;
		  txz[pos] += Ctmp * ( Dvxdz * dtdz + Dvzdx * dtdx);


    //  __syncthreads();

			zcom =  ((ix >= nabs) && (ix < NX - nabs) && (iz >= NZ - nabs)) ? 1:0;
			xcom =  (ix <= nabs) ? 0: ( (ix < NX - nabs) ? (ix - nabs) : (NX - 2*nabs));
			pmlpos = pos - ( NZ - fnabs - nabs) * (zcom + xcom);
			// update pml
			if(ix < nabs) {
				ipml = nabs - 1 - ix;
				// txx & tzz
				Vtmp = (2.0 * pml_vxx[pmlpos] + pml_para[ipml*12+5] * Dvxdx) / (2.0 + pml_para[ipml*12+1]);
				txx[pos] = txx[pos] + C[posc + 1] * dtdx * ( (pml_para[ipml*12+3] - 1) * Dvxdx - pml_para[ipml*12+3] * Vtmp * dx );
				tzz[pos] = tzz[pos] + C[posc + 2] * dtdx * ( (pml_para[ipml*12+3] - 1) * Dvxdx - pml_para[ipml*12+3] * Vtmp * dx);
				pml_vxx[pmlpos] = 2.0 * Vtmp - pml_vxx[pmlpos];
				// txz
				Vtmp = (2.0 * pml_vzx[pmlpos] + pml_para[ipml*12+4] * Dvzdx) / (2.0 + pml_para[ipml*12+0]);
				txz[pos] = txz[pos] + Ctmp * dtdx * ( (pml_para[ipml*12+2] - 1) * Dvzdx - pml_para[ipml*12+2] * Vtmp * dx );
				pml_vzx[pmlpos] = 2.0 * Vtmp - pml_vzx[pmlpos];
			}
			if(ix >= NX - nabs) {
				ipml = ix - (NX - nabs);
				// txx & tzz
				Vtmp = (2.0 * pml_vxx[pmlpos] + pml_para[ipml*12+5] * Dvxdx) / (2.0 + pml_para[ipml*12+1]);
				txx[pos] = txx[pos] + C[posc + 1] * dtdx * ( (pml_para[ipml*12+3] - 1) * Dvxdx - pml_para[ipml*12+3] * Vtmp * dx );
				tzz[pos] = tzz[pos] + C[posc + 2] * dtdx * ( (pml_para[ipml*12+3] - 1) * Dvxdx - pml_para[ipml*12+3] * Vtmp * dx);
				pml_vxx[pmlpos] = 2.0 * Vtmp - pml_vxx[pmlpos];
			}
			if(ix > NX - nabs) {
				// txz
				ipml = ix - (NX - nabs) - 1;
				Vtmp = (2.0 * pml_vzx[pmlpos] + pml_para[ipml*12+4] * Dvzdx) / (2.0 + pml_para[ipml*12+0]);
				txz[pos] = txz[pos] + Ctmp * dtdx * ( (pml_para[ipml*12+2] - 1) * Dvzdx - pml_para[ipml*12+2] * Vtmp * dx );
				pml_vzx[pmlpos] = 2.0 * Vtmp - pml_vzx[pmlpos];
			}

			if(iz < fnabs) {
				ipml = fnabs - 1 - iz;
				// txx & tzz
				Vtmp = (2.0 * pml_vzz[pmlpos] + pml_para[ipml*12+10] * Dvzdz) / (2.0 + pml_para[ipml*12+6]);
				txx[pos] = txx[pos] + C[posc + 2] * dtdz * ( (pml_para[ipml*12+8] - 1) * Dvzdz - pml_para[ipml*12+8] * Vtmp * dz );
				tzz[pos] = tzz[pos] + C[posc + 1] * dtdz * ( (pml_para[ipml*12+8] - 1) * Dvzdz - pml_para[ipml*12+8] * Vtmp * dz );
				pml_vzz[pmlpos] = 2.0 * Vtmp - pml_vzz[pmlpos];
				// txz
				Vtmp = (2.0 * pml_vxz[pmlpos] + pml_para[ipml*12+11] * Dvxdz) / (2.0 + pml_para[ipml*12+7]);
				txz[pos] = txz[pos] + Ctmp * dtdz * ( (pml_para[ipml*12+9] - 1) * Dvxdz - pml_para[ipml*12+9] * Vtmp * dz );
				pml_vxz[pmlpos] = 2.0 * Vtmp - pml_vxz[pmlpos];
			}

			if(iz >= NZ - nabs) {
				// txz
				ipml = iz - (NZ - nabs);
				Vtmp = (2.0 * pml_vxz[pmlpos] + pml_para[ipml*12+11] * Dvxdz) / (2.0 + pml_para[ipml*12+7]);
				txz[pos] = txz[pos] + Ctmp * dtdz * ( (pml_para[ipml*12+9] - 1) * Dvxdz - pml_para[ipml*12+9] * Vtmp * dz );
				pml_vxz[pmlpos] = 2.0 * Vtmp - pml_vxz[pmlpos];
			}
			if(iz > NZ - nabs) {
			// txx & tzz
				ipml = iz - (NZ - nabs) - 1;
				Vtmp = (2.0 * pml_vzz[pmlpos] + pml_para[ipml*12+10] * Dvzdz) / (2.0 + pml_para[ipml*12+6]);
				txx[pos] = txx[pos] + C[posc + 2] * dtdz * ( (pml_para[ipml*12+8] - 1) * Dvzdz - pml_para[ipml*12+8] * Vtmp * dz );
				tzz[pos] = tzz[pos] + C[posc + 1] * dtdz * ( (pml_para[ipml*12+8] - 1) * Dvzdz - pml_para[ipml*12+8] * Vtmp * dz );
				pml_vzz[pmlpos] = 2.0 * Vtmp - pml_vzz[pmlpos];
			}

    //   __syncthreads();

  		}

	return;
}
