#ifndef _MV_HOST_CU
#define _MV_HOST_CU

#include "cholmod.h"
#include "block.h"

extern "C"
void substitute_CK_host(float *L_h, size_t L_h_nz, float*bp, float *xp, size_t n);

extern "C"
void substitute_setup(float *L_h, size_t L_h_nz, float *&L_d, float *b, float *x, float *&b_x_d, size_t n);

extern "C"
void substitute_copy_back(cholmod_dense *x, double *b_x_d, size_t index);

extern "C"
void substitute_CK_free(float *L_d, double *b_x_d);

#endif
