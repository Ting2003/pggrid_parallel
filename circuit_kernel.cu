#include "circuit_kernel.h"

// kernel function, doing forward and backward substitution
// data stored in L_d: nz, col, row, and x
// data stored in b_x_d: b and x 
__global__ void substitute_CK_kernel(cudaArray *L_d, double *b_x_d){
		
}
