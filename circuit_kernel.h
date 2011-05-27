#ifndef _CIRCUIT_KERNEL_H
#define _CIRCUIT_KERNEL_H

__global__ void CK_block_kernel(float *L_d, float *b_x_d, int * L_nz_d, int *L_n_d, size_t *base_nz_d, size_t *base_n_d, int max_block_size);
#endif
