#include "triplet.h"
#include "trip_L.h"
#include <stdio.h>
#include <stdlib.h>
#include "cutil_inline.h"
#include "circuit_host.h"
#include "circuit_kernel.h"
#include "block.h"
#include "global.h"
#include <iostream>
#include <fstream>
using namespace std;

texture <float, 1, cudaReadModeElementType> L_tex;
cudaChannelFormatDesc channelDesc = 
		cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);

// kernel function, doing forward and backward substitution
// data stored in L_d: nz, col, row, and x
// data stored in b_x_d: b and x 
__global__ void substitute_CK_kernel(float *L_d, size_t L_h_nz, double *b_x_d, size_t n){
	int tid = threadIdx.x;
	// load data into shared memory
	extern __shared__ double b_x_s[];
	// two pointers to b and x
	double * b_s = b_x_s;
	double * x_s = &b_x_s[n];

	int i, j;
	int iter = n << 1 / blockIdx.x + !(n << 1 % blockIdx.x);
	for(i=0; i< iter; i++){
		int base = i * blockIdx.x;
		if((base+tid) < 2 * n)
			b_x_s[base+tid] = b_x_d[base+tid];
	}
	__syncthreads();
	
	i = 0; j = 0;
	int index_col = 0, index_row = 0;
	size_t row_p;
	// tid < WARPSIZE will do substitution
	// then all threads will copy the solution 
	// from shared memory into global memory
	if(tid < WARPSIZE){
		// doing forward substitution
		while(i < L_h_nz){
			row_p = tex1Dfetch(L_tex, i);
		
			// xj = bj / Ajj
			index_row = tex1Dfetch(L_tex, i);
			x_s[index_row] = b_s[index_row] / tex1Dfetch(L_tex, i+2);
			
			j = i+3;
			if(j >= L_h_nz) break;
			while(tex1Dfetch(L_tex, j) != tex1Dfetch(L_tex, j+1)){
				// bi = bi - Aij * xj
				index_row = tex1Dfetch(L_tex, j+1);
				index_col = tex1Dfetch(L_tex, j);
				b_s[index_row] -= tex1Dfetch(L_tex, j+2) * x_s[index_col];
				j += 3;
			}
			i = j;
		}
		
		// doing backward substitution
		i = L_h_nz - 3;
		while(i >= 0){
			row_p = tex1Dfetch(L_tex, i);
				
			// xi = bi / Aij
			x_s[row_p] = b_s[row_p] / tex1Dfetch(L_tex, i+2);

			j = i-3;
			if(j<0) break;
			
			while(tex1Dfetch(L_tex, j) != tex1Dfetch(L_tex, j+1)){
				// bi = bi - Aij * xj
				index_row = tex1Dfetch(L_tex, j+1);
				index_col = tex1Dfetch(L_tex, j);
				b_s[index_col] -= tex1Dfetch(L_tex, j+2) * x_s[index_row];
				j -= 3;
			}
			i = j;
		}
	}

	// after computing, copy back into global memory
	for(i=0; i< iter; i++){
		int base = i * blockIdx.x;
		if((base+tid) < 2 * n)
			b_x_d[base+tid] = b_x_s[base+tid];
	}
}
void substitute_CK_host(trip_L *L_h, size_t L_h_nz, cholmod_dense *b, cholmod_dense *&x){
	float *L_d = NULL;
	double *b_x_d = NULL;
	size_t n = b->nrow;
	// copy data from L to cudaArray L_d
	// bind L_d into texture memory
	substitute_setup(L_h, L_h_nz, L_d, b, x, b_x_d);
		
	dim3 dimGrid(1, 1);
	dim3 dimBlock(256, 1, 1);
	//dim3 dimBlock(WARPSIZE, 1, 1);
	// perform for- and back-ward substitution for each block
	// solution will be written from shared memory into global memory
	substitute_CK_kernel<<<dimGrid, dimBlock>>>(L_d, L_h_nz, b_x_d, n);

	// copy solution back from GPU into CPU
	// where CPU will perform the find_diff and updaterhs()
	substitute_copy_back(x, b_x_d, b->nrow);
	substitute_CK_free(L_d, b_x_d);
}

// copy data from host to device side
// 1. load sparse matrix L from host into global memory, including 4
//    sub-arrays: row index, column index, nz in each column, 
//    and value
// 2. combine L into texture memory
// 3. load dense array b and x from host into global 1d array
void substitute_setup(trip_L *L_h, size_t L_h_nz, float *L_d, cholmod_dense *b, cholmod_dense *&x, double *b_x_d){
	// count is the total bytes of all needed data in L
	// including in order of nz, col, row, and x
	size_t count_unit = sizeof(int)* 2;
	count_unit +=sizeof(double);
	
	size_t count = count_unit * L_h_nz;
	
	// allocate cudaArray and bind it with texture
	cudaMalloc((void **)&L_d, count);
	cudaMemcpy(L_d, L_h, count, cudaMemcpyHostToDevice);
	
	// bind L_d to texture memory
	cudaBindTexture(0, L_tex, L_d, count);
		
	// malloc b and x into a 1d array, copy from host to device
	count = sizeof(double)* 2 * b->nrow;
	cudaMalloc((void**)&b_x_d, count);
	size_t index = 0;
	cudaMemcpy(&b_x_d[0], b->x, b->nrow, cudaMemcpyHostToDevice);
	index += b->nrow;
	cudaMemcpy(&b_x_d[index], x->x, x->nrow, cudaMemcpyHostToDevice);
}

void substitute_copy_back(cholmod_dense *x_h, double *b_x_d, size_t index){
	size_t count = index; // both are nzmax
	cudaMemcpy(&x_h->x, &b_x_d[count], count, cudaMemcpyDeviceToHost);
}

void substitute_CK_free(float *L_d, double *b_x_d){
	cudaUnbindTexture(L_tex);
	cudaFree(b_x_d);
	cudaFree(L_d);
}
