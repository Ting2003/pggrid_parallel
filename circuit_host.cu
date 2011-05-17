#include "triplet.h"
#include <stdio.h>
#include <stdlib.h>
#include "cutil_inline.h"
#include "circuit_host.h"
#include "circuit_kernel.h"
#include "block.h"
#include <iostream>
#include <fstream>
using namespace std;

const unsigned int WARPSIZE = 32;
texture <float, 1, cudaReadModeElementType> L_tex;
cudaChannelFormatDesc channelDesc = 
		cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
	
void substitute_CK_host(cholmod_factor *L, cholmod_dense *b, cholmod_dense *&x){
	cudaArray *L_d = NULL;
	double *b_x_d = NULL;
	// copy data from L to cudaArray L_d
	// bind L_d into texture memory
	substitute_setup(L, L_d, b, x, b_x_d);
	
	dim3 dimGrid(1, 1);
	dim3 dimBlock(WARPSIZE, 1, 1);
	// perform for- and back-ward substitution for each block
	// solution will be written from shared memory into global memory
	substitute_CK_kernel<<<dimGrid, dimBlock>>>(L_d, b_x_d);

	// copy solution back from GPU into CPU
	// where CPU will perform the find_diff and updaterhs()
	substitute_copy_back(x, b_x_d, b->nzmax);
	substitute_CK_free(L_d, b_x_d);
}

// copy data from host to device side
// 1. load sparse matrix L from host into global memory, including 4
//    sub-arrays: row index, column index, nz in each column, 
//    and value
// 2. combine L into texture memory
// 3. load dense array b and x from host into global 1d array
void substitute_setup(cholmod_factor *L, cudaArray *L_d, cholmod_dense *b, cholmod_dense *&x){
	// count is the total bytes of all needed data in L
	// including in order of nz, col, row, and x
	size_t count = sizeof(size_t)*(L->n + (L->n+1)+ L->nzmax);
	count +=sizeof(double)*L->nzmax;
	
	// allocate cudaArray and bind it with texture
	cudaMallocArray(&L_d, &L_tex.channelDesc, count, 1);

	size_t index = 0;
	cudaMemcpyToArray(L_d, index, 0, L->nz, sizeof(size_t)*L->n, 
		cudaMemcpyHostToDevice);
	index += L->n;
	cudaMemcpyToArray(L_d, index, 0, L->p, sizeof(size_t)*(L->n+1), 
		cudaMemcpyHostToDevice);
	index += L->n + 1;
	cudaMemcpyToArray(L_d, index, 0, L->i, sizeof(size_t)*L->nzmax, 
		cudaMemcpyHostToDevice);
	index += L->nzmax;
	cudaMemcpyToArray(L_d, index, 0, L->x, sizeof(double)*L->nzmax, 
		cudaMemcpyHostToDevice);

	// bind L_d to texture memory
	cudaBindTextureToArray(L_tex, L_d, channelDesc);	
	
	// malloc b and x into a 1d array, copy from host to device
	count = sizeof(double)* 2 * b->nzmax;
	cudaMalloc((void**)&b_x_d, count);
	index = 0;
	cudaMemcpy(&b_x_d[0], b->x, b->nzmax, cudaMemcpyHostToDevice);
	index += b->nzmax;
	cudaMemcpy(&b_x_d[index], x->x, x->nzmax, cudaMemcpyHostToDevice);
}

void substitute_copy_back(cholmod_dense *x_h, double *b_x_d, size_t index){
	size_t count = index; // both are nzmax
	cudaMemcpy(&x_h->x[0], b_x_d, count, cudaMemcpyDeviceToHost);
}

void substitute_CK_free(cudaArray *L_d, double *b_x_d){
	cudaFree(b_x_d);
	cudaFree(L_d);
}
