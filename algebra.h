// ----------------------------------------------------------------//
// Filename : algebra.h
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// linear algebra library
// ----------------------------------------------------------------//
// - Zigang Xiao - Wed Jan 26 18:15:02 CST 2011
//   * file created

#ifndef __ALGEBRA_H__
#define __ALGEBRA_H__

#include "global.h"
#include "trip_L.h"
#include "vec.h"
#include "umfpack.h"
#include "cholmod.h"
class Vec;
class Algebra{
public:
	static void solve(const Matrix & A, const Vec & b, Vec & x);
	static void solve_CK(Matrix & A, cholmod_dense *&x, 
			cholmod_dense *b, cholmod_common *cm, size_t &peak_mem, size_t &CK_mem);
	static bool compare_row_first(const trip_L &a, const trip_L &b);
	static void factor_to_triplet(cholmod_factor *L, vector<trip_L>&L_trip);
	static void trip_to_array(vector<trip_L>&L_trip, trip_L *L_h, size_t &L_h_nz);
	static void LU_decomposition(int n, UF_long * Ap, UF_long * Ai, 
			double * Ax,
		void ** Numeric);
	static void CK_decomp(Matrix &A, cholmod_factor *&L, 
			cholmod_common *cm, size_t &peak_mem, size_t &CK_mem);
};

#endif
