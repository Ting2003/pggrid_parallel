#include "triplet.h"
#include "trip_L.h"
#include <stdio.h>
#include <stdlib.h>
#include <cutil_inline.h>
#include <cuda.h>
#include <math.h>
#include "circuit_host.h"
#include <circuit_kernel.h>
#include "block.h"
#include "global.h"
#include <iostream>
#include <fstream>
using namespace std;

texture <float, 1, cudaReadModeElementType> L_tex;
cudaChannelFormatDesc channelDesc = 
		cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);

// kernel function, doing forward and backward substitution
// data stored in L_d: row, col, and val
// data stored in b_x_d: b and x 
/*__global__ void substitute_CK_kernel(float *L_d, size_t L_h_nz, float *b_x_d, size_t n){
	int tid = threadIdx.x;
	// load data into shared memory
	extern __shared__ float b_x_s[];

	int i, j;
	int iter = n / blockDim.x ;
	if((n % blockDim.x)!=0) iter += 1;
	for(i=0; i< iter; i++){
		int base = i * blockDim.x;
		if((base+tid) < n)
			b_x_s[base+tid] = b_x_d[base+tid];
	}
	__syncthreads();
		
	i = 0; j = 0;
	int index_col = 0, index_row = 0;
	size_t row_p;
	// tid < WARPSIZE will do substitution
	// then all threads will copy the solution 
	// from shared memory into global memory
	if(tid < HALF_WARP){
		// doing forward substitution
		while(i < 3*L_h_nz){
			row_p = tex1Dfetch(L_tex, i);
		
			// xj = bj / Ajj
			index_row = tex1Dfetch(L_tex, i);
			b_x_s[index_row] /= tex1Dfetch(L_tex, i+2);
			
			j = i+3;
			if(j >= 3 * L_h_nz) break;
			//while(L_d[j] != L_d[j+1]){
			while(tex1Dfetch(L_tex, j) != tex1Dfetch(L_tex, j+1)){
				// bi = bi - Aij * xj
				index_row = tex1Dfetch(L_tex, j);
				index_col = tex1Dfetch(L_tex, j+1);
				b_x_s[index_row] -= tex1Dfetch(L_tex, j+2) * b_x_s[index_col];
				j += 3;
			}
			i = j;
		}
		
		// doing backward substitution
		i = 3 * L_h_nz - 3;
		while(i >= 0){
			row_p = tex1Dfetch(L_tex, i);
				
			// xi = bi / Aij
			b_x_s[row_p] /= tex1Dfetch(L_tex, i+2);

			j = i-3;
			if(j<0) break;
			//while(L_d[j]!=L_d[j+1]){
			while(tex1Dfetch(L_tex, j) != tex1Dfetch(L_tex, j+1)){
				// bi = bi - Aij * xj
				index_row = tex1Dfetch(L_tex, j);
				index_col = tex1Dfetch(L_tex, j+1);
				b_x_s[index_col] -= tex1Dfetch(L_tex, j+2) * b_x_s[index_row];
				j -= 3;
			}
			i = j;
		}
	}
		
	// after computing, copy back into global memory
	for(i=0; i< iter; i++){
		int base = i * blockDim.x;
		if((base+tid) < n)
			//b_x_d[base+tid] = tex1Dfetch(L_tex, tid);
			b_x_d[base+tid] = b_x_s[base+tid];
	}
}*/

void substitute_CK_host(float *L_h, size_t L_h_nz, float *bp, float *xp, size_t n){
	float *L_d ;	
	float *b_x_d ;
	// dump cudaMalloc, as the first call wll cost about 1s
	// which is much larger than usual 1e-6s
	cudaMalloc((void**)&L_d, sizeof(float));
	unsigned int timer;
	float cudaTime;
	CUT_SAFE_CALL(cutCreateTimer(&timer));
	CUT_SAFE_CALL(cutStartTimer(timer));

	// copy data from L to cudaArray L_d
	// bind L_d into texture memory	
	substitute_setup(L_h, L_h_nz, L_d, bp, xp, b_x_d, n);

	dim3 dimGrid(1, 1);
	dim3 dimBlock(256, 1, 1);
	int sharedMemSize =sizeof(float) *n;
	clog<<"shared mem size: "<<sharedMemSize<<endl;
	// perform for- and back-ward substitution for each block
	// solution will be written from shared memory into global memory
	unsigned int timer_compute;
	float cudaTime_compute;
	CUT_SAFE_CALL(cutCreateTimer(&timer_compute));
	CUT_SAFE_CALL(cutStartTimer(timer_compute));
	
	//substitute_CK_kernel<<<dimGrid, dimBlock, sharedMemSize>>>(L_d, L_h_nz, b_x_d, n);
	
	cutilCheckMsg("Kernel execution failed.");
	CUT_SAFE_CALL(cutStopTimer(timer_compute));
	cudaTime_compute = cutGetTimerValue(timer_compute);
	clog<<"kernel time: "<<cudaTime_compute/1000<<" (s) "<<endl;
	CUT_SAFE_CALL(cutDeleteTimer(timer_compute));	

	// copy solution back from GPU into CPU
	// where CPU will perform the find_diff and updaterhs()
	//cutilSafeCall(cudaMemcpy(bp, b_x_d, sizeof(float)*n, cudaMemcpyDeviceToHost));
	cutilSafeCall(cudaMemcpy(xp, b_x_d, sizeof(float)*n, cudaMemcpyDeviceToHost));
	
	cutilSafeCall(cudaUnbindTexture(L_tex));
	cutilSafeCall(cudaFree(L_d));
	cutilSafeCall(cudaFree(b_x_d));
	
	CUT_SAFE_CALL(cutStopTimer(timer));
	cudaTime = cutGetTimerValue(timer);
	clog<<"gpu time: "<<cudaTime/1000<<" (s) "<<endl;
	CUT_SAFE_CALL(cutDeleteTimer(timer));
}


// copy data from host to device side
// 1. load sparse matrix L from host into global memory, including 4
//    sub-arrays: row index, column index, nz in each column, 
//    and value
// 2. combine L into texture memory
// 3. load dense array b and x from host into global 1d array
void substitute_setup(float *L_h, size_t L_h_nz, float *&L_d, float *bp, float *xp, float *&b_x_d, size_t n){
	size_t count = sizeof(float) * 3 * L_h_nz;		
	// allocate cudaArray and bind it with texture
	cudaMalloc((void**)&L_d, count);		
	clog<<"Texture memory size: "<<count<<endl;
	cutilSafeCall(cudaMemcpy(L_d, L_h, count, cudaMemcpyHostToDevice));	
	cutilSafeCall(cudaBindTexture(0, L_tex, L_d, channelDesc, count));
		
	// malloc b and x into a 1d array, copy from host to device
	count = sizeof(float) * n;
	cutilSafeCall(cudaMalloc((void**)&b_x_d, count));	
	cutilSafeCall(cudaMemcpy(b_x_d, bp, count, cudaMemcpyHostToDevice));
}

// block version of substitution
void block_CK_host(BlockInfo &block_info){
	float *L_d ;
	float *b_x_d ;
		
	// dump cudaMalloc, as the first call wll cost about 1s
	// which is much larger than usual 1e-6s
	cudaMalloc((void**)&L_d, sizeof(float));
	
	unsigned int timer;
	float cudaTime;
	CUT_SAFE_CALL(cutCreateTimer(&timer));
	CUT_SAFE_CALL(cutStartTimer(timer));
	
	// copy data from L to cudaArray L_d
	// bind L_d into texture memory
	// copy L_n_nz_h into constant memory in 1D array
	size_t total_n = 0;
	size_t total_nz = 0;
	int *L_nz_d = NULL; size_t *base_nz_d = NULL;
	int *L_n_d = NULL; size_t *base_n_d = NULL;
	substitute_block_setup(block_info, L_d, b_x_d, L_nz_d, 
		L_n_d, base_nz_d, base_n_d, total_n, total_nz);
	
	// find max block size in block_info
	int max_block_size = 0;
	for(size_t i=0;i<block_info.size();i++)
		if(block_info[i].count > max_block_size)
			max_block_size = block_info[i].count;

	dim3 dimGrid(block_info.X_BLOCKS, block_info.Y_BLOCKS);
	dim3 dimBlock(256, 1, 1);
	
	// calculate the shared memory size
	// note: shared memory is in unit size, max
	int sharedMemSize =sizeof(float) *max_block_size;
	clog<<"shared mem size: "<<sharedMemSize<<endl;
	// perform for- and back-ward substitution for each block
	// solution will be written from shared memory into global memory
	unsigned int timer_compute;
	float cudaTime_compute;
	CUT_SAFE_CALL(cutCreateTimer(&timer_compute));
	CUT_SAFE_CALL(cutStartTimer(timer_compute));
	
	CK_block_kernel<<<dimGrid, dimBlock, sharedMemSize>>>
		(L_d, b_x_d, L_nz_d, L_n_d, base_nz_d, base_n_d, max_block_size);
	
	cutilCheckMsg("Kernel execution failed.");
	CUT_SAFE_CALL(cutStopTimer(timer_compute));
	cudaTime_compute = cutGetTimerValue(timer_compute);
	clog<<"kernel time: "<<cudaTime_compute/1000<<" (s) "<<endl;
	CUT_SAFE_CALL(cutDeleteTimer(timer_compute));	

	// copy solution back from GPU into CPU
	// where CPU will perform the find_diff and updaterhs()
	size_t base = 0;
	for(size_t i=0;i<block_info.size();i++){
		//clog<<"block index: "<<i<<endl;
		cudaMemcpy(block_info[i].xp_f, 
		&b_x_d[base], sizeof(float)*block_info[i].count, 
		cudaMemcpyDeviceToHost);
		base  += block_info[i].count;
	}
	//cutilSafeCall(cudaMemcpy(xp, b_x_d, sizeof(float)*n, cudaMemcpyDeviceToHost));
	
	cutilSafeCall(cudaUnbindTexture(L_tex));
	cudaFree(L_d);
	cudaFree(b_x_d);
	
	CUT_SAFE_CALL(cutStopTimer(timer));
	cudaTime = cutGetTimerValue(timer);
	clog<<"gpu time: "<<cudaTime/1000<<" (s) "<<endl;
	CUT_SAFE_CALL(cutDeleteTimer(timer));
}

// copy data from host to device side
// 1. load sparse matrix L from host into global memory, 
//    including 3 sub-arrays: row index, column index, and value
// 2. combine L into texture memory
// 3. load dense array b and x from host into global 1d array
// 4. copy L_n_nz indo constant memory
void substitute_block_setup(BlockInfo &block_info, float *&L_d, 
	float *&b_x_d, int *&L_nz_d, int *&L_n_d, 
	size_t *&base_nz_d, size_t *&base_n_d, size_t &total_n, 
	size_t &total_nz){
	// L_nz_h stores nz, L_n_h stores n for each block
	int *L_nz_h; size_t *base_nz_h;
	int *L_n_h; size_t *base_n_h;

	L_nz_h = new int [block_info.size()];
	L_n_h = new int [block_info.size()];
	base_nz_h = new size_t [block_info.size()];
	base_n_h = new size_t [block_info.size()];

	for(size_t i=0;i<block_info.size();i++){
		L_n_h[i] = block_info[i].count;
		L_nz_h[i] = block_info[i].L_h_nz;
		if(i==0){
			base_nz_h[i] = 0;
			base_n_h[i] = 0;
		}
		else{
			base_nz_h[i] = base_nz_h[i-1] +
				       3*block_info[i-1].L_h_nz;
			base_n_h[i] = base_n_h[i-1]+
				       block_info[i-1].count;
		} 
	}
	total_n = base_n_h[block_info.size()-1] + 
		  block_info[block_info.size()-1].count;
	total_nz = base_nz_h[block_info.size()-1] + 
		  block_info[block_info.size()-1].L_h_nz;

	// Allocate and copy L_nz_d, as well as L_n_d
	size_t count = sizeof(int)*block_info.size();
	cutilSafeCall(cudaMalloc((void**)&L_nz_d, count));
	cutilSafeCall(cudaMalloc((void**)&L_n_d, count));
	cutilSafeCall(cudaMemcpy(L_nz_d, L_nz_h, count, 
			cudaMemcpyHostToDevice));
	cutilSafeCall(cudaMemcpy(L_n_d, L_n_h, count, 
			cudaMemcpyHostToDevice));
	
	// Allocate and copy base_nz_d, as well as base_n_d
	count = sizeof(size_t) *block_info.size();
	cutilSafeCall(cudaMalloc((void**)&base_nz_d, count));
	cutilSafeCall(cudaMalloc((void**)&base_n_d, count));
	cutilSafeCall(cudaMemcpy(base_nz_d, base_nz_h, count, 
			cudaMemcpyHostToDevice));
	cutilSafeCall(cudaMemcpy(base_n_d, base_n_h, count, 
			cudaMemcpyHostToDevice));
	
	// free host memory
	free(L_nz_h); free(L_n_h);
	free(base_nz_h); free(base_n_h);
		
	// allocate cudaArray and bind it with texture
	count = sizeof(float) * 3 * total_nz;
	cutilSafeCall(cudaMalloc((void**)&L_d, count));		
	clog<<"Texture memory size: "<<count<<endl;
	size_t base = 0;
	for(size_t i=0;i<block_info.size();i++){
		cutilSafeCall(cudaMemcpy(&L_d[base], block_info[i].L_h, 
		sizeof(float)*3*block_info[i].L_h_nz, cudaMemcpyHostToDevice));
		base += 3*block_info[i].L_h_nz;
	}
	cutilSafeCall(cudaBindTexture(0, L_tex, L_d, channelDesc, count));
		
	// malloc b and x into a 1d array, copy from host to device
	count = sizeof(float) * total_n;
	cutilSafeCall(cudaMalloc((void**)&b_x_d, count));
	base = 0;
	for(size_t i=0;i<block_info.size();i++){
		//clog<<"block id: "<<block_info[i].bid<<endl;
		//for(size_t j=0;j<block_info[i].count;j++)
			//clog<<j<<" "<<block_info[i].bnewp_f[j]<<endl;
		cutilSafeCall(cudaMemcpy(&b_x_d[base], block_info[i].bnewp_f, 
		sizeof(float)*block_info[i].count, cudaMemcpyHostToDevice));
		base += block_info[i].count;
	}
}
