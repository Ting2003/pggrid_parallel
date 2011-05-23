#include "circuit_kernel.h"
#include "global.h"

// kernel function, doing forward and backward substitution
// data stored in L_d: nz, col, row, and x
// data stored in b_x_d: b and x 
__global__ void substitute_CK_kernel(float *L_d, size_t L_h_nz, float *b_x_d, size_t n){
	int tid = threadIdx.x;
	// load data into shared memory
	extern __shared__ float b_x_s[];
	// two pointers to b and x
	float * b_s = b_x_s;
	float * x_s = &b_x_s[n];

	int i, j;
	int iter = (n << 1) / blockDim.x ;
	if(((n << 1) % blockDim.x)!=0) iter += 1;
	for(i=0; i< iter; i++){
		int base = i * blockDim.x;
		if((base+tid) < 2 * n){
			b_x_s[base+tid] = b_x_d[base+tid]+1;
			b_x_s[base+tid] += 1;
		}
			//b_x_s[base+tid] = b_x_d[base+tid];
	}
	__syncthreads();
	
	/*i = 0; j = 0;
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
	*/	
	// after computing, copy back into global memory
	for(i=0; i< iter; i++){
		int base = i * blockDim.x;
		if((base+tid) < 2 * n)
			b_x_d[base+tid] = b_x_s[base+tid];
	}
}
