#ifndef _ARITHMETIC_GRID2_PRECONDITIONER_BLOCK_REDBLACK_H_
#define _ARITHMETIC_GRID2_PRECONDITIONER_BLOCK_REDBLACK_H_

#include "../../../Basis/Basis/Environments.h"
#include "Grid2_PoissonMatrix.h"

namespace arithmetic {

template<typename T>
class Grid2_Preconditioner_BlockRedBlackICCG : public Grid2_Preconditioner<T> {
public : 
	Grid2_Preconditioner_BlockRedBlackICCG() : matrix(NULL) {}
	Grid2_Preconditioner_BlockRedBlackICCG(grid_index_type nX, grid_index_type nY, grid_index_type nbX, grid_index_type nbY) : matrix(NULL) { resize(nX,nY); makeBlocks(nbX, nbY); }

	void resize(grid_index_type, grid_index_type);

	void makeBlocks(index_type nbX, index_type nbY);

	void buildPreconditioner(const Grid2_PoissonMatrix<T>& mat);
	void applyPreconditioner(Grid2<T,T>&, const Grid2<T,T>&) const;

public : 
	const Grid2_PoissonMatrix<T>* matrix;

	Grid2<T,T>				preconditioner;
	mutable Grid2<T,T>		m;

	struct Block {
		index_type begin[2], end[2];
		Block(index_type b[2], index_type e[2]) { begin[0] = b[0]; begin[1] = b[1]; end[0] = e[0]; end[1] = e[1]; }		
	};

private : 
	std::vector<Block>		red;
	std::vector<Block>		black;
};

template<typename T> 
void Grid2_Preconditioner_BlockRedBlackICCG<T>::resize(grid_index_type nX, grid_index_type nY) {
	matrix = NULL;
	preconditioner.init(nX, nY, 0, 0);
	m.init(nX, nY, 0, 0);
	red.clear(); black.clear();
}

template<typename T> 
void Grid2_Preconditioner_BlockRedBlackICCG<T>::makeBlocks(index_type nbX, index_type nbY) {
	index_type nX = preconditioner.nX, nY = preconditioner.nY;
	if(nX == 0 || nY == 0) 
		return;

	red.clear(); black.clear();

	index_type nx = (nX-2) / nbX;	
	index_type ny = (nY-2) / nbY;

	if(nx * nbX != nX-2) nx++;
	if(ny * nbY != nY-2) ny++;

	index_type begin[2], end[2];
	for(index_type i=0;i<nbX;++i) {
		for(index_type j=0;j<nbY;++j) {
			begin[0] = i*nx+1;	end[0] = MIN(begin[0]+nx, nX-1);
			begin[1] = j*ny+1;	end[1] = MIN(begin[1]+ny, nY-1);

			if((i+j) % 2 == 0)	red.push_back(Block(begin, end));
			else				black.push_back(Block(begin, end));
		}
	}
}

template<typename T>
void Grid2_Preconditioner_BlockRedBlackICCG<T>::buildPreconditioner(const Grid2_PoissonMatrix<T>& mat) {
	matrix = &mat;

	preconditioner.zero();

	BuildPreconditioner_RedBlock_Grid2<T> func_red(mat.isBoundary, mat.poisson, preconditioner, red);
	func_red(basis::blocked_range<index_type>(0, index_type(red.size())));

	BuildPreconditioner_BlackBlock_Grid2<T> func_black(mat.isBoundary, mat.poisson, preconditioner, black);
	func_black(basis::blocked_range<index_type>(0, index_type(black.size())));
}

template<typename T>
void Grid2_Preconditioner_BlockRedBlackICCG<T>::applyPreconditioner(Grid2<T,T> &y, const Grid2<T,T> &x) const {
	m.zero(); y.zero();

	ForwardSubstitution_RedBlock_Grid2<T> func_fw_red(matrix->isBoundary, matrix->poisson, preconditioner, red, x, m);
	func_fw_red(basis::blocked_range<index_type>(0, index_type(red.size())));

	ForwardSubstitution_BlackBlock_Grid2<T> func_fw_black(matrix->isBoundary, matrix->poisson, preconditioner, black, x, m);
	func_fw_black(basis::blocked_range<index_type>(0, index_type(black.size())));

	BackwardSubstitution_BlackBlock_Grid2<T> func_bw_black(matrix->isBoundary, matrix->poisson, preconditioner, black, m, y);
	func_bw_black(basis::blocked_range<index_type>(0, index_type(black.size())));

	BackwardSubstitution_RedBlock_Grid2<T> func_bw_red(matrix->isBoundary, matrix->poisson, preconditioner, red, m, y);
	func_bw_red(basis::blocked_range<index_type>(0, index_type(red.size())));
}

};		// arithmetic

#include "Grid2_Preconditioner_BlockRedBlackICCG_RedFunctors.h"
#include "Grid2_Preconditioner_BlockRedBlackICCG_BlackFunctors.h"

#endif