#ifndef _MV_HOST_CU
#define _MV_HOST_CU

#include "cholmod.h"
#include "block.h"

extern "C"
void substitute_CK_host(float *L_h, size_t L_h_nz, float*bp, float *xp, size_t n);

extern "C"
void substitute_setup(float *L_h, size_t L_h_nz, float *&L_d, float *b, float *x, float *&b_x_d, size_t n);

// block version function
extern "C" 
void block_CK_host(BlockInfo &block_info);

extern "C"
void substitute_block_setup(BlockInfo &block_info, float *&L_d, float *&b_x_d, int *&L_nz_d, int *&L_n_d, size_t *&base_nz_d, size_t *&base_n_d, size_t &toal_n, size_t &total_nz);


#endif
