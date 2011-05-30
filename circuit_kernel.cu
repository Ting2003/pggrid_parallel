#include "circuit_kernel.h"
#include "global.h"

// block version kernel function 
// doing forward and backward substitution
// data stored in L_d: row, col, val for each block
// data stored in b_x_d: b for each block
// L_nz_d, L_n_d: nz and n for each block
// base_nz_d, base_n_d: nz base and n base for each block
__global__ void CK_block_kernel(float *L_d, float *b_x_d, 
		int * L_nz_d, int *L_n_d, size_t *base_nz_d, 
		size_t *base_n_d, int max_block_size){
	int tid = threadIdx.x;
	// load 1 block data into shared memory
	extern __shared__ float b_x_s[];

	int block_id = blockIdx.y * gridDim.x + blockIdx.x;
	
	long i, j;
	int iter = L_n_d[block_id] / blockDim.x ;
	if(((L_n_d[block_id] % blockDim.x)!=0)) iter += 1;
	size_t block_base_n = base_n_d[block_id];
	size_t block_base_nz = base_nz_d[block_id];
	for(i=0; i< iter; i++){
		int thread_base = i * blockDim.x;
		if((thread_base+tid) < L_n_d[block_id])
			b_x_s[thread_base+tid] = b_x_d[thread_base+tid+block_base_n];
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
		while(i < 3*L_nz_d[block_id]){
			row_p =L_d[block_base_nz+i];
		
			// xj = bj / Ajj
			index_row = L_d[block_base_nz+i];
			b_x_s[index_row] /= L_d[block_base_nz+i+2];
			
			j = i+3;
			if(j >= 3 * L_nz_d[block_id]-2) break;
			while(L_d[block_base_nz+j] != L_d[block_base_nz+j+1]){
				// bi = bi - Aij * xj
				index_row = L_d[block_base_nz+j];
				index_col = L_d[block_base_nz+j+1];
				b_x_s[index_row] -= L_d[block_base_nz+j+2]*b_x_s[index_col];
				j += 3;
			}
			i = j;
		}
			
		// doing backward substitution
		i = 3 * L_nz_d[block_id] - 3;
		while(i >= 0){
			row_p = L_d[block_base_nz+i];
				
			// xi = bi / Aij
			b_x_s[row_p] /= L_d[block_base_nz+i+2];

			j = i-3;
			if(j<0) break;
			while(L_d[block_base_nz+j]!=L_d[block_base_nz+j+1]){
				// bi = bi - Aij * xj
				index_row = L_d[block_base_nz+j];
				index_col = L_d[block_base_nz+j+1];
				b_x_s[index_col] -= L_d[block_base_nz+j+2] * b_x_s[index_row];
				j -= 3;
			}
			i = j;
		}
	}
	// after computing, copy back into global memory
	for(i=0; i< iter; i++){
		int thread_base = i * blockDim.x;
		if((thread_base+tid) < L_n_d[block_id])
			b_x_d[thread_base+tid+block_base_n] = b_x_s[tid+thread_base];
	}
}
