#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"


float* init_matrix(int sizenum){
	int i;
	float *U;
	if( ! (U = (float *)malloc(sizenum*sizeof(float))) )
      printf("init failed !\n");
	for (i = 0; i < sizenum; i++) U[i] = 0;
	return U;
}

int* init_matrix_int(int sizenum){
	int i;
	int *U;
	U = (int *)malloc((sizenum*sizeof(int)));
	for (i = 0; i < sizenum; i++) U[i] = 0;
	return U;
}

float Ricker(float time, float freq0, float amp0 ){
	float wave, a;
	time = time - ( 1.5/freq0 ); // delay
	//time = time - ( 1.5/freq0 ); // delay
  a = freq0 * sqrt(2.) * pi;
  // if (time > tsrcf) wave = 0.;
  wave = (1. - a*a*time*time)*exp(-0.5*a*a*time*time);// Ricker
  //wave = a*a*time*(a*a*time*time-3.0)*exp(-0.5*a*a*time*time);  // Ricker derivative
  wave = wave * amp0;
  return(wave);
}


void set_pml(int nabs, float dx, float dz, float dt, float Vmax, float freq0, float* pml_para){

	int i, id;
	float lnR, thk_pml, beta0, alpha0, d0, xL;
	float d0factor = 3.0;
	float PPW0=10.0;
	float p_power=2.0;

	// natural logarithm of the theoretical reflection coefficient R
	lnR = log(10) * (-3.0 - (log10( nabs + 0.0 ) - 1.0) / log10(2.0));
	// thickness of PML

  for(id = 0; id < 2; id++) {
    if(id == 0) dx = dx;
    if(id == 1) dx = dz;
    //printf("dx %f\n", dx);
	  thk_pml = nabs * dx;
	  d0 = -(p_power + 1.0) * Vmax * lnR / ( 2.0 * thk_pml);
	  d0 = d0 * d0factor;
	  beta0 = Vmax / (0.5 * PPW0 * dx * freq0);
	  if(beta0 < 1.0)  beta0 = 1.0;
	  alpha0 = pi * freq0;

	  for(i = 0; i < nabs; i++) {
	  	// i=0: interior interface; i=nabs-1: exterior boundry;

	  	// define damping profile at grid points
	  	xL=(i+1)*dx/thk_pml;
	  	pml_para[i*12+id*6+4]=d0*pow(xL,p_power);
	  	pml_para[i*12+id*6+2]=1.0+(beta0-1.0)*pow(xL,p_power);
	  	pml_para[i*12+id*6+0]=alpha0*(1.0-xL);

	  	// define damping profile at half grid points
	  	//xL=(i+1)*dx/thk_pml;
	  	xL=(i+0.5)*dx/thk_pml;
	  	pml_para[i*12+id*6+5]=d0*pow(xL,p_power);
	  	pml_para[i*12+id*6+3]=1.0+(beta0-1.0)*pow(xL,p_power);
	    pml_para[i*12+id*6+1]=alpha0*(1.0-xL);

	    if(pml_para[i*12+id*6+0]<0.0)  pml_para[i*12+id*6+0]=0.0;
	  	if(pml_para[i*12+id*6+1]<0.0)  pml_para[i*12+id*6+1]=0.0;

	  	// beta <-- 1/beta
	  	pml_para[i*12+id*6+2]=1.0/pml_para[i*12+id*6+2];
	  	pml_para[i*12+id*6+3]=1.0/pml_para[i*12+id*6+3];

	  	// d <-- d/beta
	  	pml_para[i*12+id*6+4]=pml_para[i*12+id*6+4]*pml_para[i*12+id*6+2];
	  	pml_para[i*12+id*6+5]=pml_para[i*12+id*6+5]*pml_para[i*12+id*6+3];

	  	// alpha <-- alpha + d/beta
	  	pml_para[i*12+id*6+0]=pml_para[i*12+id*6+0]+pml_para[i*12+id*6+4];
	  	pml_para[i*12+id*6+1]=pml_para[i*12+id*6+1]+pml_para[i*12+id*6+5];

	  	// multiply alpha+d/beta by dt
	  	pml_para[i*12+id*6+0]=dt*pml_para[i*12+id*6+0];
	  	pml_para[i*12+id*6+1]=dt*pml_para[i*12+id*6+1];
	  	// multiply d/beta by dt/dx
	  	pml_para[i*12+id*6+4]=dt/dx*pml_para[i*12+id*6+4];
	  	pml_para[i*12+id*6+5]=dt/dx*pml_para[i*12+id*6+5];
    }
  }

  return;
}

void get_C(FILE* mdfid, int NK, int NZ, int NX, int nabs, int fnabs, int readin_model,\
    float* C, float& Vmax, float& Vmin){

  //float Vmax = 0.;
  //float Vmin = 1e10;
  float rho = 2.200;
  float vptmp = 2000;
  float vstmp = 1200;
  int ix, iz, ik;
  float s1, s2, s3, tmpf;
  //for( ix = nabs; ix < NX - nabs; ix++ )
    for( iz = fnabs; iz < NZ - nabs; iz++ )
      for( ix = nabs; ix < NX - nabs; ix++ ) {
      if(readin_model){
   		  fread(&tmpf,sizeof(float),1,mdfid); s1=tmpf;//vp
   		  fread(&tmpf,sizeof(float),1,mdfid); s2=tmpf;//vs
   		  fread(&tmpf,sizeof(float),1,mdfid); s3=tmpf*0.001;//rho
        C[ix * NK * NZ + iz * NK + 0] = s3;//rho;//s3;
        C[ix * NK * NZ + iz * NK + 1] = s3*s1*s1;//rho*vptmp*vptmp*dt/dx;//lam2mu;
        C[ix * NK * NZ + iz * NK + 2] = s3*s1*s1 - 2*s3*s2*s2;//;//lam;
        C[ix * NK * NZ + iz * NK + 3] = s3*s2*s2;//rho*vstmp*vstmp*dt/dx;//mu;
        if (s1 > Vmax) Vmax = s1;
        if (s2 < Vmin) Vmin = s2;
      }
      else{
        C[ix * NK * NZ + iz * NK + 0] = rho;//s3;
        C[ix * NK * NZ + iz * NK + 1] = rho*vptmp*vptmp;//lam2mu;
        C[ix * NK * NZ + iz * NK + 2] = rho*vptmp*vptmp - 2*rho*vstmp*vstmp;//lam;
        C[ix * NK * NZ + iz * NK + 3] = rho*vstmp*vstmp;//mu;
        if (s1 > Vmax) Vmax = vptmp;//s1;
        if (s2 < Vmin) Vmin = vstmp;// s2;

      }
   	}
	fclose(mdfid);
	//printf("%f\n", C[nabs*NZ*4 + fnabs*4 + 0]);
    //boundary
  for(ix = 0; ix < nabs; ix++)
  	for(iz = fnabs; iz < NZ - nabs; iz++)
  		for(ik = 0; ik < 4; ik++ )
  			C[ix * NK * NZ + iz * NK + ik] = C[nabs * NK * NZ + iz * NK + ik];
  for(ix = NX - nabs; ix < NX; ix++ )
   	for(iz = fnabs; iz < NZ - nabs; iz++ )
   		for(ik = 0; ik < 4; ik++ )
   			C[ix * NK * NZ + iz * NK + ik] = C[(NX - nabs - 1) * NK * NZ + iz * NK + ik];
  for(ix = 0; ix < NX; ix++ )
   	for(iz = 0; iz < fnabs; iz++ )
   		for(ik = 0; ik < 4; ik++ )
   			C[ix * NK * NZ + iz * NK + ik] = C[ix * NK * NZ + fnabs * NK + ik];
  for(ix = 0; ix < NX; ix++ )
   	for(iz = NZ - nabs; iz < NZ; iz++ )
   		for(ik = 0; ik < 4; ik++ )
   			C[ix * NK * NZ + iz * NK + ik] = C[ix * NK * NZ + (NZ - nabs - 1) * NK + ik];

}


void readin_coor(FILE* srfid, int N, int* loc, float dx, float dz){

  int ii;
  float tmpx, tmpz;
  for (ii=0; ii < N; ii++){
    fscanf(srfid, "%f %f",&tmpx, &tmpz);
    loc[ii*2+0] = int(tmpx / dx );
    loc[ii*2+1] = int(tmpz / dz );
  }

  fclose(srfid);
  return;
}

void mymemorycp(float* buf_array, float* tar, int start, int len){

	int ii;
	for (ii = 0; ii < len; ii++){
		buf_array[ii] =  tar[start + ii];
	}
	return;
}
