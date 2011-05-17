#ifndef _MV_HOST_CU
#define _MV_HOST_CU

#include "cholmod.h"
#include "block.h"

extern "C"
void substitute_CK_host(cholmod_factor *L, cudaArray *L_d, cholmod_dense *b, cholmod_dense *&x);

extern "C"
void substitute_setup(cholmod_factor *L, cudaArray *L_d, cholmod_dense *b, cholmod_dense *&x);
#endif
