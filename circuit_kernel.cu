#include "circuit_kernel.h"
#include "global.h"
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
	int index_col = 0;
	size_t row_p, row_q;
	float sum = 0;
	// tid < WARPSIZE will do substitution
	// then all threads will copy the solution 
	// from shared memory into global memory
	if(tid < WARPSIZE){
		// doing forward substitution
		while(i < L_h_nz){
			j = i+3;
			if(j >= L_h_nz) break;
			row_p = tex1Dfecth(L_tex, i);
			col_p = tex1Dfetch(L_tex, i+1);
			if(row_p == row_q){
				// sum = Aij*xj
				index_col = tex1Dfetch(L_tex, i+1);
				sum += x_s[index_col] * tex1Dfetch(L_tex, i+2); 	
			}
			else if(row_p ==  col_p){
				// xi = (bi - sum) / Aii
				index_col = tex1Dfecth(L_tex, i);
				x_s[index_col] = (b_s[index_col] - sum) / tex1Dfecth(L_tex, i+2);
				sum = 0;
			}
			i += 3;
		}
		// last element for last row
		// xi = (bi - sum) / Aii
		i = L_h_nz - 3;
		index_col = tex1Dfecth(L_tex, i);
		x_s[index_col] = (b_s[index_col] - sum) / tex1Dfecth(L_tex, i+2);
		sum = 0;

		// doing backward substitution
		i = L_h_nz - 3;
		while(i >= 0){
			j = i-3;
			if(j<0) break;
			row_p = tex1Dfetch(L_tex, i);
			row_q = tex1Dfetch(L_tex, j);
			if(row_p == row_q){
				// sum+ = Aij * xj
				index_col = tex1Dfetch(L_tex, i+1);
				sum += x_s[index_col] * tex1Dfetch(L_tex, i+2);		
			}
			else{
				// xi = (bi - sum) / Aii
				index_col = tex1Dfetch(L_tex, i);
				x_s[index_col] = (b_s[index_col] - sum) / tex1Dfetch(L_tex, i+2);
				sum = 0;
			}
			i -= 3;
		}
		// last element for first row
		// xi = (bi - sum) / Aii
		i = 0;
		index_col = tex1Dfecth(L_tex, i);
		x_s[index_col] = (b_s[index_col] - sum) / tex1Dfecth(L_tex, i+2);
		sum = 0;
	}

}
