#ifndef _ARITHMETIC_GRID3_PRECONDITIONER_BLOCK_REDBLACK_H_
#define _ARITHMETIC_GRID3_PRECONDITIONER_BLOCK_REDBLACK_H_

#include "Grid3_PoissonMatrix.h"

namespace arithmetic {

template<typename T>
class Grid3_Preconditioner_BlockRedBlackICCG : public Grid3_Preconditioner<T> {
public : 
	Grid3_Preconditioner_BlockRedBlackICCG() : matrix(NULL) {}
	Grid3_Preconditioner_BlockRedBlackICCG(index_type nX, index_type nY, index_type nZ, index_type nbX, index_type nbY, index_type nbZ) : matrix(NULL) { resize(nX,nY,nZ); makeBlocks(nbX,nbY,nbZ); }

	void resize(index_type, index_type, index_type);

	void makeBlocks(index_type nbX, index_type nbY, index_type nbZ);

	void buildPreconditioner(const Grid3_PoissonMatrix<T>& mat);
	void applyPreconditioner(Grid3<T,T>&, const Grid3<T,T>&) const;

public : 
	const Grid3_PoissonMatrix<T>* matrix;

	Grid3<T,T>				preconditioner;
	mutable Grid3<T,T>		m;

	struct Block {
		index_type begin[3], end[3];
		Block(index_type _b[3], index_type _e[3]) { for(index_type i=0;i<3;++i) { begin[i] = _b[i]; end[i] = _e[i]; } }		
	};

private : 
	std::vector<Block>		red;
	std::vector<Block>		black;
};

template<typename T> 
void Grid3_Preconditioner_BlockRedBlackICCG<T>::resize(index_type nX, index_type nY, index_type nZ) {
	matrix = NULL;
	preconditioner.init(nX, nY, nZ, 0, 0, 0);
	m.init(nX, nY, nZ, 0, 0, 0);
	red.clear(); black.clear();
}

template<typename T> 
void Grid3_Preconditioner_BlockRedBlackICCG<T>::makeBlocks(index_type nbX, index_type nbY, index_type nbZ) {
	index_type nX = preconditioner.nX, nY = preconditioner.nY, nZ = preconditioner.nZ;
	if(nX == 0 || nY == 0 || nZ == 0) 
		return;

	red.clear(); black.clear();

	index_type nx = (preconditioner.nX-2) / nbX;	
	index_type ny = (preconditioner.nY-2) / nbY;
	index_type nz = (preconditioner.nZ-2) / nbZ;

	if(nx * nbX != preconditioner.nX-2) nx++;
	if(ny * nbY != preconditioner.nY-2) ny++;
	if(nz * nbZ != preconditioner.nZ-2) nz++;

	index_type begin[3], end[3];
	for(index_type i=0;i<nbX;++i) {
		for(index_type j=0;j<nbY;++j) {
			for(index_type k=0;k<nbZ;++k) {
				begin[0] = i*nx+1;	end[0] = MIN(begin[0]+nx, nX-1);
				begin[1] = j*ny+1;	end[1] = MIN(begin[1]+ny, nY-1);
				begin[2] = k*nz+1;	end[2] = MIN(begin[2]+nz, nZ-1);

				if((i+j+k) % 2 == 0)	red.push_back(Block(begin, end));
				else					black.push_back(Block(begin, end));
			}			
		}
	}
}

template<typename T>
void Grid3_Preconditioner_BlockRedBlackICCG<T>::buildPreconditioner(const Grid3_PoissonMatrix<T>& mat) {
	matrix = &mat;

	preconditioner.zero();

	BuildPreconditioner_RedBlock_Grid3<T> func_red(mat.isBoundary, mat.poisson, preconditioner, red);
	func_red(basis::blocked_range<index_type>(0, index_type(red.size())));

	BuildPreconditioner_BlackBlock_Grid3<T> func_black(mat.isBoundary, mat.poisson, preconditioner, black);
	func_black(basis::blocked_range<index_type>(0, index_type(black.size())));
}

template<typename T>
void Grid3_Preconditioner_BlockRedBlackICCG<T>::applyPreconditioner(Grid3<T,T> &y, const Grid3<T,T> &x) const {
	m.zero(); y.zero();

	ForwardSubstitution_RedBlock_Grid3<T> func_fw_red(matrix->isBoundary, matrix->poisson, preconditioner, red, x, m);
	func_fw_red(basis::blocked_range<index_type>(0, index_type(red.size())));

	ForwardSubstitution_BlackBlock_Grid3<T> func_fw_black(matrix->isBoundary, matrix->poisson, preconditioner, black, x, m);
	func_fw_black(basis::blocked_range<index_type>(0, index_type(black.size())));

	BackwardSubstitution_BlackBlock_Grid3<T> func_bw_black(matrix->isBoundary, matrix->poisson, preconditioner, black, m, y);
	func_bw_black(basis::blocked_range<index_type>(0, index_type(black.size())));

	BackwardSubstitution_RedBlock_Grid3<T> func_bw_red(matrix->isBoundary, matrix->poisson, preconditioner, red, m, y);
	func_bw_red(basis::blocked_range<index_type>(0, index_type(red.size())));
}

};			// arithmetic

#include "Grid3_Preconditioner_BlockRedBlackICCG_RedFunctors.h"
#include "Grid3_Preconditioner_BlockRedBlackICCG_BlackFunctors.h"

#endif
