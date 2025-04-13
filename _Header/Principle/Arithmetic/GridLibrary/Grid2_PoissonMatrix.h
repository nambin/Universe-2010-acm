#ifndef _ARITHMETIC_GRID2_POISSON_MATRIX_H_
#define _ARITHMETIC_GRID2_POISSON_MATRIX_H_

#include "../../../Basis/Basis/Environments.h"
#include "../../../Basis/ParallelProgramming/ParallelComputable.h"

namespace arithmetic {

template<typename T>
class Grid2_PoissonMatrix : public ParallelComputable {
public : 
	typedef typename Grid2<T,T>::index_type index_type;

	Grid2_PoissonMatrix() {}
	Grid2_PoissonMatrix(index_type nX, index_type nY, const T h, const Vector2<T>& offset) { resize(nX, nY, h, offset); }

	index_type nX() const { return isBoundary.nX; }
	index_type nY() const { return isBoundary.nY; }
	index_type size() const { return isBoundary.size(); }

	void resize(index_type, index_type, const T h, const Vector2<T>& offset);
	template<typename R> void resize(const Grid2<R,T>& grid) { resize(grid.nX, grid.nY, grid.dx, grid.offset); }

	template<typename VECTOR>
	void mul(Grid2<T,T>&, const VECTOR&) const;

	const bool importBinary(std::ifstream&);
	const bool exportBinary(std::ofstream&) const;

public : 
	Grid2<bool,T>		isBoundary;
	Grid2x3<T>			poisson;
};

template<typename T>
const bool Grid2_PoissonMatrix<T>::importBinary(std::ifstream& ifs) {
	isBoundary.importBinary(ifs);
	poisson.importBinary(ifs);

	return true;
}

template<typename T>
const bool Grid2_PoissonMatrix<T>::exportBinary(std::ofstream& ofs) const {
	isBoundary.exportBinary(ofs);
	poisson.exportBinary(ofs);

	return true;
}

template<typename T>
void Grid2_PoissonMatrix<T>::resize(index_type nX, index_type nY, const T h, const Vector2<T>& offset) {
	isBoundary.init(nX, nY, h, h, offset);
	poisson.init(nX, nY);

	isBoundary.set(false);
}

template<typename T> template<typename GRID2>
void Grid2_PoissonMatrix<T>::mul(Grid2<T,T>& y, const GRID2& x) const {
	y.zero();

#ifdef USE_OPENMP
	int my_end = int(y.nX-1);
	int i;

#pragma omp parallel for
#else
	index_type my_end = y.nX-1;
	index_type i;
#endif
	for( i = 1 ; i < my_end; ++i ) {
		for( index_type j = 1; j < y.nY-1; ++j ) {
			if( !isBoundary(i,j) )	{
				y(i,j) =  poisson(i,j,0)   * x(i,j)
						+ poisson(i-1,j,1) * x(i-1,j)
						+ poisson(i,j,1)   * x(i+1,j)
						+ poisson(i,j-1,2) * x(i,j-1)
						+ poisson(i,j,2)   * x(i,j+1);	
			}				
		}
	}

}

template<typename T>
class Grid2_Preconditioner {
public : 
	typedef typename Grid2<T,T>::index_type index_type;
	
	Grid2_Preconditioner() {}

	virtual void buildPreconditioner(const Grid2_PoissonMatrix<T>& mat) = 0;
	virtual void applyPreconditioner(Grid2<T,T>&, const Grid2<T,T>&) const = 0;
};

template<typename T>
class Grid2_Preconditioner_Identity : public Grid2_Preconditioner<T> {
public : 
	Grid2_Preconditioner_Identity() {}

	void buildPreconditioner(const Grid2_PoissonMatrix<T>& mat) {}
	void applyPreconditioner(Grid2<T,T>& x, const Grid2<T,T>& y) const { x = y; }
};

};		// arithmetic

#endif