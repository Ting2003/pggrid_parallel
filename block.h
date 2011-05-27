// ----------------------------------------------------------------//
// Filename : block.h
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// declaration of block class
// ----------------------------------------------------------------//

#ifndef __BLOCK_H__
#define __BLOCK_H__
#include <fstream>
#include "triplet.h"
#include "global.h"
#include "vec.h"
#include "net.h"
#include "util.h"
#include "umfpack.h"
#include "cholmod.h"

using namespace std;

class Block{
	typedef vector<Net *> NetPtrVector;
public:
	Block(size_t count=0);
	~Block();
	void free_block_cholmod(cholmod_common *cm);
	void LU_decomposition();
	void CK_decomp(Matrix & A, cholmod_common *cm, size_t &peak_mem, size_t &CK_mem);
	void solve_CK(cholmod_common *cm); // solve with cholesky decomp
	void solve_CK_setup(cholmod_common *cm);

	// allocate space for the matrix and vectors accoding to count
	void allocate_resource(cholmod_common*cm);

	DIRECTION relative_to(const Block & block) const;

	void update_x();

	void update_rhs();

	bool inside_bbox(long x, long y) const{
		return (x>=lx && x<=ux && y>=ly && y<=uy);
	}

	cholmod_factor * L;
	// variable for parallel usage
	float *L_h;
	size_t L_h_nz; // #nz in L

	NetPtrVector boundary_netlist;
	
	// vector b
	cholmod_dense * b_ck, *b_new_ck;
	// pointer to b_ck, b_new_ck, and x_ck;
	double *bp, *bnewp, *xp;

	// float pointer for parallel computation
	float *bnewp_f, *xp_f;
	// solution
	cholmod_dense *x_ck;

	// its id in the block_info
	size_t bx, by;
	size_t bid;

	// number of *representative* nodes in this block
	// equal to matrix size and b size
	size_t count;

	Node ** nodes;

	// geometric information of this block
	double lx, ly, ux, uy;
};

class BlockInfo{
public:
	BlockInfo():X_BLOCKS(1), Y_BLOCKS(1),
	len_per_block_x(0), len_per_block_y(0),
	len_ovr_x(0), len_ovr_y(0){};

	// update the right hand side according to the boundary netlist and
	// value of x
	void update_rhs(size_t block_id, vector<DIRECTION> & dirs);
	void update_rhs(size_t block_id); // just a shorthand for updating four dirs

	Block * get_block_neighbor(const Block & b, DIRECTION dir);

	void resize(size_t n){blocks.resize(n);}

	Block & operator [] (size_t index){ return blocks[index]; }

	size_t size() const{return blocks.size();}

	void update_block_geometry();

	// set the x, y length per block
	void set_len_per_block(size_t x_min, size_t x_max,
			       size_t y_min, size_t y_max,
			       double overlap_ratio);

	size_t X_BLOCKS; // # of blocks along x axis
	size_t Y_BLOCKS; // # of blocks along y axis

	double len_per_block_x; // length per block along x direction
	double len_per_block_y; // length per block along y direction
	double len_ovr_x; // overlap length along x direction
	double len_ovr_y; // overlap length along y direction
private:
	vector< Block > blocks;
};

#endif
