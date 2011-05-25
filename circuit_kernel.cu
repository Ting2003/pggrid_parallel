#include "circuit_kernel.h"
#include "global.h"

//extern texture<float, 1> L_tex;
// kernel function, doing forward and backward substitution
// data stored in L_d: nz, col, row, and x
// data stored in b_x_d: b and x 
__global__ void substitute_CK_kernel(float *L_d, size_t L_h_nz, float *b_x_d, size_t n){
	int tid = threadIdx.x;
	// load data into shared memory
	extern __shared__ float b_x_s[];

	int i, j;
	int iter = n / blockDim.x ;
	if(((n % blockDim.x)!=0)) iter += 1;
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
			row_p =L_d[i];
		
			// xj = bj / Ajj
			index_row = L_d[i];
			b_x_s[index_row] /= L_d[i+2];
			
			j = i+3;
			if(j >= 3 * L_h_nz-2) break;
			while(L_d[j] != L_d[j+1]){
				// bi = bi - Aij * xj
				index_row = L_d[j];
				index_col = L_d[j+1];
				b_x_s[index_row] -= L_d[j+2]*b_x_s[index_col];
				j += 3;
			}
			i = j;
		}
			
		// doing backward substitution
		i = 3 * L_h_nz - 3;
		while(i >= 0){
			row_p = L_d[i];
				
			// xi = bi / Aij
			b_x_s[row_p] /= L_d[i+2];

			j = i-3;
			if(j<0) break;
			while(L_d[j]!=L_d[j+1]){
				// bi = bi - Aij * xj
				index_row = L_d[j];
				index_col = L_d[j+1];
				b_x_s[index_col] -= L_d[j+2] * b_x_s[index_row];
				j -= 3;
			}
			i = j;
		}
	}
	// after computing, copy back into global memory
	for(i=0; i< iter; i++){
		int base = i * blockDim.x;
		if((base+tid) < n)
			b_x_d[base+tid] = b_x_s[base+tid];
	}
}
