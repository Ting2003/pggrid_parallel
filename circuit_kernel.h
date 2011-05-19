#ifndef _CIRCUIT_KERNEL_H
#define _CIRCUIT_KERNEL_H
__global__ void substitute_CK_kernel(float *L_d, size_t L_h_nz, double *b_x_d, size_t n);
#endif
