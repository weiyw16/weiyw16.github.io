/**************************************
***		Forward modeling			***
*** 	Author: Yanwen Wei          ***
*** 	Date: 2019-9-23             ***
***		Email: wei_yanwen@163.com   ***
***************************************/


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "common.h"
#include "kernel.cuh"


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

void copy_float(float* ori, float* dev_ori, int msize){

  cudaError_t err = cudaSuccess;
  err = cudaMemcpy(dev_ori, ori, msize*sizeof(float), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) {
    printf("stderr, %s\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }

}
void copy_int(int* ori, int* dev_ori, int msize){

  cudaError_t err = cudaSuccess;
  err = cudaMemcpy(dev_ori, ori, msize*sizeof(int), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) {
    printf("stderr, %s\n", cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }
}






int main(int argc, char** argv){

// parameters
	time_t t0;
	clock_t start, stop;
  int readin_model, save_snap, readin_source, readin_receiver;
  int snap_step;
	int is, ir, it, ix, iz, nabs, fnabs, pos;
	int nx, nz, NX, NZ, NK, NS, NR, NT;
	int sz, sx0, sdx, rz0, rdz, rx;
	int *src_loc, *rec_loc, *this_src_loc;
	float dx, dz, dt, freq0, amp0;
	float Vmin, Vmax, wave;
	float *record_vz, *record_vx, *div, *curl, *snap_vz, *snap_vx;
	float *C, *txx, *tzz, *txz, *vx, *vz;// readin C
	float *pml_vxx, *pml_vzz, *pml_vzx, *pml_vxz;//NX*NZ - nx*nnz
	float *pml_xtxx, *pml_xtxz, *pml_ztxz, *pml_ztzz;//NX*NZ - nx*nnz
  float *pml_para;// nabs*12
	char logfile[25], paramfile[25], modelfile[25];
  char sourcefile[25], receiverfile[25], snapf_vz[126], snapf_vx[126];
	char outfile_vz[80], outfile_vx[80], outfile_div[80], outfile_curl[80];
	FILE *fid, *pafid, *mdfid, *oufid, *srfid;


	// default value of the parameters
	// ****************************************************************
  readin_model = 0; save_snap = 0; 
  readin_source = 0; readin_receiver = 0;
	sprintf(logfile,"log_m%s_s%s.txt", argv[1], argv[2]);
//	sprintf(paramfile, "ParamInput.txt");
  sprintf(paramfile, "ParamInput_%s.txt", argv[2]);
	//sprintf(modelfile, "modelfile.bin");
	// sprintf(outfile_vz, "vz.bin");
	// sprintf(outfile_vx, "vx.bin");
	// sprintf(outfile_div, "div.bin");
	// sprintf(outfile_curl, "curl.bin");
	//sprintf(sourcefile, "sourcefile.txt");
	//sprintf(receiverfile, "receiverfile.txt");
	nx = 300; nz = 400;
	nabs = 40; fnabs = 40;
	NX = nx + 2 * nabs; NZ = nz + fnabs + nabs;
	NS = 1; NR = 0; NT = 100; NK = 4;
	dx = 5; dz = 5; dt = 0.001; freq0 = 20; amp0 = 1;
  snap_step = int(0.1 / dt);
	// ****************************************************************

	time(&t0);//local time


// prepare the log file
	// sprintf(logfile,"log.txt");
	if( !(fid=fopen(logfile,"w")) ) 
		{ printf("log file not opened"); return(1); }
	fprintf(fid, "***************************************************\n");
	fprintf(fid, "\tBegin to Produce Synthetic Seismic Data \n");
	fprintf(fid, "\tstart time: %s",ctime(&t0));
	fprintf(fid, "***************************************************\n\n");

// readin parameters
	if( !(pafid=fopen(paramfile,"r")) ) {
		fprintf(fid,"Error opening param file %s\n",paramfile);}
	else{
		fprintf(fid,"Start reading parameters in file %s\n", paramfile);
		fscanf(pafid,"%i %i %f %f", &nx, &nz, &dx, &dz);
		fscanf(pafid,"%i %i %i %i", &sz, &sx0, &sdx, &NS);
		fscanf(pafid,"%i %i %i %i", &rx, &rz0, &rdz, &NR);
		fscanf(pafid,"%f %f %f %i", &freq0, &amp0, &dt, &NT);
		fscanf(pafid,"%i %i", &nabs, &fnabs);
		fscanf(pafid,"%i", &readin_model);
		fscanf(pafid,"%i %i", &save_snap, &snap_step);
		fscanf(pafid,"%i %i", &readin_source, &readin_receiver);
	}
	fclose(pafid);
// update parameters and initial matrix
	fprintf(fid, "***\t\tThe running parameters\t\t***\n");
  fprintf(fid, "nx %i, nz %i dx %f dz %f\n", nx, nz, dx, dz);
  fprintf(fid, "sz %i, sx0 %i sdx %i NS %i\n", sz, sx0, sdx, NS);
  fprintf(fid, "rx %i, rz0 %i rdz %i NR %i\n", rx, rz0, rdz, NR);
  fprintf(fid, "freq0 %f, amp0 %f dt %f NT %i\n", freq0, amp0, dt, NT);
  fprintf(fid, "nabs %i, fnabs %i\n", nabs, fnabs);
  fprintf(fid, "save_snap %i, snap_step %i\n", save_snap, snap_step);
  fprintf(fid, "readin_source %i, readin_receiver %i\n", readin_source, readin_receiver);
	// C = init_matrix( NX * NZ * NK);


// init all the global matrix
	NX = nx + 2 * nabs; NZ = nz + fnabs + nabs;
	C = init_matrix( NX * NZ * NK);
	vx = init_matrix( NX * NZ );
	vz = init_matrix( NX * NZ );
	txx = init_matrix( NX * NZ );
	tzz = init_matrix( NX * NZ );
	txz = init_matrix( NX * NZ );

	pml_vxx = init_matrix( NX * NZ - nx * nz);
	pml_vzz = init_matrix( NX * NZ - nx * nz);
	pml_vxz = init_matrix( NX * NZ - nx * nz);
	pml_vzx = init_matrix( NX * NZ - nx * nz);
	pml_xtxx = init_matrix( NX * NZ - nx * nz);
	pml_xtxz = init_matrix( NX * NZ - nx * nz);
	pml_ztxz = init_matrix( NX * NZ - nx * nz);
	pml_ztzz = init_matrix( NX * NZ - nx * nz);

  pml_para = init_matrix( nabs*12 );// 6 * 2

	src_loc = init_matrix_int( 2 * NS);
	rec_loc = init_matrix_int( 2 * NR);
	this_src_loc = init_matrix_int( 2 );
	record_vz = init_matrix( NR * NT );
	record_vx = init_matrix( NR * NT );
	div = init_matrix( NR * NT );
	curl = init_matrix( NR * NT );
	snap_vz = init_matrix( nx * nz );
	snap_vx = init_matrix( nx * nz );


  // set the source and receiver locations
  if (! readin_source){
	  for(is = 0; is < NS; is++){
	  	src_loc[2 * is + 0] = sx0 + sdx * is;
	  	src_loc[2 * is + 1] = sz;
	  }
  }
  else{
    sprintf(sourcefile, "sourcefile.txt");
    if( !(srfid=fopen(sourcefile,"r")) ) {
    fprintf(fid,"Error opening param file %s\n", sourcefile);
    	return(1); }
    else{
     fprintf(fid,"\n\nStart reading sourcefile %s\n", sourcefile);
     readin_coor(srfid, NS, src_loc, dx, dz);
    } 
  }
  if (! readin_receiver){
  	for(ir = 0; ir < NR; ir++){
  		rec_loc[2 * ir + 0] = rx;
  		rec_loc[2 * ir + 1] = rz0 + rdz * ir;
  	}
  }
  else{
    sprintf(receiverfile, "receiverfile.txt");
    if( !(srfid=fopen(receiverfile,"r")) ) {
    fprintf(fid,"Error opening param file %s\n", receiverfile);
    	return(1); }
    else{
     fprintf(fid,"\n\nStart reading receiverfile %s\n", receiverfile);
     readin_coor(srfid, NR, rec_loc, dx, dz);
    }
  }

  fprintf(fid,"\n\n souce position \n\n");
  for(ix=0; ix<NS; ix++){fprintf(fid, "%i %i\n", src_loc[ix*2+0], src_loc[ix*2+1]);}
  fprintf(fid,"\n\n receiver position \n\n");
  for(ix=0; ix<NR; ix++){fprintf(fid, "%i %i\n", rec_loc[ix*2+0], rec_loc[ix*2+1]);}

 // readin model and initial C
 sprintf(modelfile, "modelfile_%s.bin", argv[1]);
 if( !(mdfid=fopen(modelfile,"rb")) ) {
 fprintf(fid,"Error opening param file %s\n", modelfile);
 	return(1); }
 else{
  fprintf(fid,"\n\nStart reading modelfile %s\n", modelfile);
  Vmax = 0.; Vmin = 1e10;
  get_C (mdfid, NK, NZ, NX, nabs, fnabs, readin_model, C, Vmax, Vmin);
 }

// check
  fprintf(fid, "***CHECKPOINT***\n");
  fprintf(fid, "Vmax, %f\n", Vmax);
  fprintf(fid, "Vmix, %f\n", Vmin);
  fprintf(fid, "VZ, %f\n", vz[100]);

  set_pml(nabs, dx, dz, dt, Vmax, freq0, pml_para);



  cudaSetDevice(atoi(argv[3]));

  int pml_size = NZ*NX - nz*nx;
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
  float *dev_pml_para = get_deviceMem_float(pml_para, nabs*12);
  int *dev_rec_loc = get_deviceMem_int(rec_loc, NR);
  int *dev_src_loc = get_deviceMem_int(src_loc, 2);
  int block_size = 32;
  dim3 Threads(block_size, block_size, 1);
  //dim3 Grids((NX + block_size - 1) / block_size, (NZ + block_size - 1) / block_size, 1);
  dim3 Grids((NX - 1) / block_size + 1 , (NZ - 1) / block_size + 1, 1);

// kernel 
  time(&t0);//local time
  fprintf(fid, "\n\n***************************************************\n");
  fprintf(fid, "\t loops start at time %s", ctime(&t0));
  fprintf(fid, "***************************************************");

  for (is = 0; is < NS; is++){

  	this_src_loc[0] = src_loc[is*2];
  	this_src_loc[1] = src_loc[is*2 + 1];
    copy_int(this_src_loc, dev_src_loc, 2);
    //int *dev_src_loc = get_deviceMem_int(this_src_loc, 2);
//    checkCudaErrors( cudaMemcpy(dev_src_loc, this_src_loc, 2*sizeof(int), cudaMemcpyHostToDevice) );
    cudaDeviceSynchronize();
  	start = clock();

  	for (it = 0; it < NT; it++){

  		wave = Ricker(it*dt, freq0, amp0);
      wave = - wave * dt / (dx * dx);
  //		printf("it %d, wave %f\n", it, wave);
      
  		///update_stress(vx, vz, txx, tzz, txz,\
			  C, pml_vxx, pml_vzz, pml_vxz, pml_vzx, pml_para,\
			  NZ, NX, NK, nabs, fnabs,\
			  wave, this_src_loc, dx, dz, dt);
      kernel_update_stress <<< Grids, Threads >>> ( dev_vx,  dev_vz,  dev_txx,  dev_tzz,  dev_txz,\
        dev_C, dev_pml_vxx,  dev_pml_vzz,  dev_pml_vxz,  dev_pml_vzx, dev_pml_para,\
	      NZ,  NX,  NK,  nabs,  fnabs,\
	      wave, dev_src_loc, dx,  dz,  dt);
      cudaDeviceSynchronize();
     // copyback_float(txx, dev_txx, NX*NZ);
     // copyback_float(tzz, dev_tzz, NX*NZ);
     // copyback_float(txz, dev_txz, NX*NZ);
     // cudaDeviceSynchronize();
	    //update_velocity(vx, vz, txx, tzz, txz,\
			  C, pml_vxx, pml_vzz, pml_vxz, pml_vzx,\
			  pml_xtxx, pml_xtxz, pml_ztxz, pml_ztzz, pml_para,\
			  NZ, NX, NK, NR, nabs, fnabs,\
			  rec_loc, record_vz, record_vx, div, curl,\
        wave, this_src_loc, it, dx, dz, dt);
      //copy_float(vx, dev_vx, NX*NZ);
      //copy_float(vz, dev_vz, NX*NZ);

	    //gpu_update_velocity(vx, vz, txx, tzz, txz,\
			  C, pml_vxx, pml_vzz, pml_vxz, pml_vzx,\
			  pml_xtxx, pml_xtxz, pml_ztxz, pml_ztzz, pml_para,\
			  NZ, NX, NK, NR, nabs, fnabs,\
			  rec_loc, record_vz, record_vx, div, curl,\
        wave, this_src_loc, it, dx, dz, dt, NX*NZ-nz*nx);
      kernel_update_velocity <<< Grids, Threads >>> ( dev_vx,  dev_vz,  dev_txx,  dev_tzz,  dev_txz,\
       dev_C,	dev_pml_xtxx,  dev_pml_xtxz,  dev_pml_ztxz,  dev_pml_ztzz,  dev_pml_para,\
       dev_rec_loc,  dev_src_loc,  wave,\
	     NZ,  NX,  NK,  NR,  nabs,  fnabs,\
	     dx,  dz,  dt, it);
/*  
     //kernel_update <<< Grids, Threads >>> ( dev_vx,  dev_vz,  dev_txx,  dev_tzz,  dev_txz,\
       dev_C, dev_pml_vxx,  dev_pml_vzz,  dev_pml_vxz,  dev_pml_vzx,\
	     dev_pml_xtxx,  dev_pml_xtxz,  dev_pml_ztxz,  dev_pml_ztzz,  dev_pml_para,\
       dev_rec_loc,  dev_src_loc,  wave,\
	     NZ,  NX,  NK,  NR,  nabs,  fnabs,\
	     dx,  dz,  dt);
       */
     cudaDeviceSynchronize();
     copyback_float(vx, dev_vx, NX*NZ);
     copyback_float(vz, dev_vz, NX*NZ);

   
	  for(ir = 0; ir < NR; ir++){
	  	pos = (rec_loc[ir * 2 + 0] + nabs) * NZ + rec_loc[ir * 2 + 1] + fnabs;
//  		record_vx[it * NR + ir] = vx[pos];
//  		record_vz[it * NR + ir] = vz[pos];
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

    if (save_snap && it % snap_step == 0){
      for (ix=0; ix < nx; ix++)
        for (iz=0; iz < nz; iz++){
          snap_vz[ix*nz + iz] = vz[(ix+nabs)*NZ + iz + fnabs];
          snap_vx[ix*nz + iz] = vx[(ix+nabs)*NZ + iz + fnabs];
      }
  	  sprintf(snapf_vz, "snap_output/model_%s_shot_%s_it_%d_vz.bin", argv[1], argv[2], it);
  	  if( !(oufid=fopen(snapf_vz,"wb")) ) {
		  fprintf(fid,"Error opening output file %s\n", snapf_vz);
      	return(1); }
  	  else fwrite(snap_vz, sizeof(float), nz*nx, oufid);
  	  fclose(oufid);
  	  sprintf(snapf_vx, "snap_output/model_%s_shot_%s_it_%d_vx.bin", argv[1], argv[2], it);
  	  if( !(oufid=fopen(snapf_vx,"wb")) ) {
		  fprintf(fid,"Error opening output file %s\n", snapf_vx);
      	return(1); }
  	  else fwrite(snap_vx, sizeof(float), nz*nx, oufid);
  	  fclose(oufid);
      }
  
  

       
    }// end of time iteration

  	stop = clock();
  	fprintf(fid, "\nFinish shot %d, time consuming %f sec.", is, (double)(stop - start)/CLOCKS_PER_SEC);




  	sprintf(outfile_vz, "output/model_%s_shot_%s_vz.bin", argv[1], argv[2]);
  	sprintf(outfile_vx, "output/model_%s_shot_%s_vx.bin", argv[1], argv[2]);
  	sprintf(outfile_div, "output/model_%s_shot_%s_div.bin", argv[1], argv[2]);
  	sprintf(outfile_curl, "output/model_%s_shot_%s_curl.bin", argv[1], argv[2]);

  	if( !(oufid=fopen(outfile_vz,"wb")) ) {
		fprintf(fid,"Error opening output file %s\n", outfile_vz);
    	return(1); }
  	else fwrite(record_vz, sizeof(float), NT*NR, oufid);
  	fclose(oufid);
  
  	if( !(oufid=fopen(outfile_vx,"wb")) ) {
  		fprintf(fid,"Error opening output file %s\n", outfile_vx);
      	return(1); }
  	else fwrite(record_vx, sizeof(float), NT*NR, oufid);
  	fclose(oufid);
  
  	if( !(oufid=fopen(outfile_div,"wb")) ) {
  		fprintf(fid,"Error opening output file %s\n", outfile_div);
      	return(1); }
  	else fwrite(div, sizeof(float), NT*NR, oufid);
  	fclose(oufid);
  
  	if( !(oufid=fopen(outfile_curl,"wb")) ) {
  		fprintf(fid,"Error opening output file %s\n", outfile_curl);
      	return(1); }
  	else fwrite(curl, sizeof(float), NT*NR, oufid);
  	fclose(oufid);

  }


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



// close
	free(vx);free(vz);free(txx);free(tzz);free(txz);
	free(pml_vxx);free(pml_vxz);free(pml_vzz);free(pml_vzx);
	free(pml_xtxx);free(pml_ztzz);free(pml_xtxz);free(pml_ztxz);
  free(pml_para);
	free(C);free(rec_loc);free(src_loc);free(this_src_loc);
	free(record_vz);free(record_vx);free(div);free(curl);
  free(snap_vz); free(snap_vx);
	time(&t0);//local time
    fprintf(fid, "\n\n***************************************************\n");
    fprintf(fid, "\t done at time %s", ctime(&t0));
    fprintf(fid, "***************************************************");
	fclose(fid);
	
	return 0;
}
