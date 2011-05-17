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

texture <float, 1, cudaReadModeElementType> L_tex;
cudaChannelFormatDesc channelDesc = 
		cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
	
void substitute_CK_host(cholmod_factor *L, cholmod_dense *b, cholmod_dense *&x){
	cudaArray *L_d = NULL;
	// copy data from L to host L_d
	// bind L_d into texture memory
	substitute_setup(L, L_d, b, x);
	
	/*// malloc and copy data from host to device
	Matrix A_d[block_info.size()];
	decomp_setup(A, A_d, block_info);
	//configuration part goes here
	size_t x = block_info.X_BLOCKS;
	size_t y = block_info.Y_BLOCKS;
	dim3 dimGrid(x, y);
	dim3 dimBlock(WARPSIZE,1,1);
	ck_decomp_kernel<<<dimGrid, dimBlock>>>();
	decomp_copy_back(block_info);
	// copy data back to host
	//decomp_free();
	*/
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
	double *b_x_d = NULL;
	count = sizeof(double)* 2 * b->nzmax;
	cudaMalloc((void**)&b_x_d, count);
	index = 0;
	cudaMemcpy(&b_x_d[0], b->x, b->nzmax, cudaMemcpyHostToDevice);
	index += b->nzmax;
	cudaMemcpy(&b_x_d[index], x->x, x->nzmax, cudaMemcpyHostToDevice);
}
