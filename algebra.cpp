// ----------------------------------------------------------------//
// Filename : algebra.cpp
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// linear algebra library
// ----------------------------------------------------------------//
// - Zigang Xiao - Wed Jan 26 18:15:02 CST 2011
//   * algebra::solve

#include "cholmod.h"
#include <algorithm>
#include <cassert>
#include <ctype.h>
#include "global.h"
#include "trip_L.h"
#include "util.h"
#include "vec.h"
#include "umfpack.h"
#include "algebra.h"

// order: col > row > val
bool Algebra::compare_row_first(const trip_L &a, const trip_L &b){
	if(a.row == b.row){
		return (a.col < b.col);
	}
	return (a.row < b.row);
}

// transform factor matrix from column to triplet
// output is column-wise triplet expression of L
void Algebra::factor_to_triplet(cholmod_factor *L, vector<trip_L> &L_trip){
	int *L_nz, *L_p, *L_i;
	double *L_x;
	L_nz = static_cast<int *> (L->nz);
	L_p = static_cast<int *> (L->p);
	L_i = static_cast<int *> (L->i);
	L_x = static_cast<double *> (L->x);
	size_t n = L->n;
	size_t base = 0;
	trip_L temp;
	for(size_t i=0; i< n; i++){
		for(int j=0; j< L_nz[i]; j++){
			temp.row = L_i[base+j];
			temp.col = i;
			temp.val = L_x[base+j];
			L_trip.push_back(temp);	
		}
		base += L_nz[i];
	}
	for(size_t i=0;i<L_trip.size();i++)
		clog<<"L_trip in column: "<<L_trip[i].row<<
		" "<<L_trip[i].col<<" "<<L_trip[i].val<<endl;
	//sort(L_trip.begin(), L_trip.end(), compare_row_first);
	//free(L_nz); free (L_p); free(L_i); free (L_x);
}

void Algebra::trip_to_array(vector<trip_L>&L_trip, float *&L_h, size_t &L_h_nz){
	L_h = new float[3 * L_trip.size()];
	for(size_t i=0; i<L_trip.size();i++){
		L_h[3*i] = L_trip[i].row;
		L_h[3*i+1] = L_trip[i].col;
		L_h[3*i+2] = L_trip[i].val;
	}
	// L_h_nz is the length of L_h
	L_h_nz = L_trip.size();
	L_trip.clear();
}

// deliver the address of x
void Algebra::solve_CK(Matrix & A, cholmod_dense *&x, cholmod_dense *b, cholmod_common *cm, size_t &peak_mem, size_t &CK_mem, float *bp, float *xp){
	cholmod_factor *L;
	cm->final_ll = true; //stay in LL' format
	CK_decomp(A, L, cm, peak_mem, CK_mem);
	// then solve
	// x = cholmod_solve(CHOLMOD_A, L, b, cm);
	// hello();
	vector<trip_L> L_trip;
	// L_h is the memory used for host memory, in array format
	float * L_h = NULL;
	// record the length of L_h
	size_t L_h_nz = 0;
	factor_to_triplet(L, L_trip);
	
	trip_to_array(L_trip, L_h, L_h_nz);
	for(size_t i=0;i<b->nrow;i++)
		clog<<"i, bp, xp: "<<i<<" "<<bp[i]<<" "<<xp[i]<<endl;
	substitute_CK_host(L_h, L_h_nz, bp, xp, b->nrow);
	cholmod_print_dense(b, "b", cm);
	cholmod_print_factor( L, "L",cm);
	cholmod_free_factor(&L, cm);	
}

// doing cholesky decomposition
void Algebra::CK_decomp(Matrix &A, cholmod_factor *&L, cholmod_common *cm, size_t &peak_mem, size_t & CK_mem){
	// doing factorization first
	cholmod_triplet * T;
	size_t n_row = A.get_row();
	size_t n_col = A.get_row();
	size_t nnz = A.size();
	int *Ti;
	int *Tj;
	double *Tx;
	int stype = -1;// lower triangular storage
	T = cholmod_allocate_triplet(n_row, n_col, nnz, stype, 
			CHOLMOD_REAL, cm);
	Ti = static_cast<int *>(T->i);
	Tj = static_cast<int *>(T->j);
	Tx = static_cast<double *>(T->x);
	// copy data into T
	for(size_t k=0;k<nnz;k++){
		Ti[k] = A.Ti[k];
		Tj[k] = A.Tj[k];
		Tx[k] = A.Tx[k];
	}
	T->nnz = nnz;
	A.Ti.clear();
	A.Tj.clear();
	A.Tx.clear();
	cholmod_sparse * A_cholmod;
	A_cholmod = cholmod_triplet_to_sparse(T, nnz, cm);

	// free the triplet pointer
	cholmod_free_triplet(&T, cm);

	cm->supernodal = -1;// always do simplictical LL'	
	//cm->final_super = false;
	//cm->final_ll = true;
	//cholmod_print_common("common", cm);
	L = cholmod_analyze(A_cholmod, cm);
	cholmod_factorize(A_cholmod, L, cm);
	if(peak_mem < cm->memory_usage)
		peak_mem = cm->memory_usage;
	CK_mem += cm->lnz;
	cholmod_free_sparse(&A_cholmod, cm);
}
