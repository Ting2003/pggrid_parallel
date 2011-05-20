#ifndef _MV_HOST_CU
#define _MV_HOST_CU

#include "cholmod.h"
#include "block.h"

extern "C"
void substitute_CK_host(trip_L *L_h, size_t L_h_nz, cholmod_dense *b, cholmod_dense *&x);

extern "C"
void substitute_setup(trip_L *L_h, size_t L_h_nz, float *L_d, cholmod_dense *b, cholmod_dense *&x, double *b_x_d);

extern "C"
void substitute_copy_back(cholmod_dense *x, double *b_x_d, size_t index);

extern "C"
void substitute_CK_free(float *L_d, double *b_x_d);

#endif
