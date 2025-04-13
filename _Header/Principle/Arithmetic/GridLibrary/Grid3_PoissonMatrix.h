#ifndef _ARITHMETIC_GRID3_POISSON_MATRIX_H_
#define _ARITHMETIC_GRID3_POISSON_MATRIX_H_

#include "../../../Basis/Basis/Environments.h"
#include "../../../Basis/ParallelProgramming/ParallelComputable.h"
#include "GridLibrary.h"

namespace arithmetic {

template<typename T>
class Grid3_PoissonMatrix : public ParallelComputable {
public : 
	typedef typename Grid3<T,T>::index_type index_type;

	Grid3_PoissonMatrix() {}
	Grid3_PoissonMatrix(index_type nX, index_type nY, index_type nZ, const T h, const Vector3<T>& offset) { resize(nX, nY, nZ, h, offset); }

	index_type nX() const { return isBoundary.nX; }
	index_type nY() const { return isBoundary.nY; }
	index_type nZ() const { return isBoundary.nZ; }
	index_type size() const { return isBoundary.size(); }

	void resize(index_type, index_type, index_type, const T h, const Vector3<T>& offset);
	template<typename R> void resize(const Grid3<R,T>& grid) { resize(grid.nX, grid.nY, grid.nZ, grid.dx, grid.offset); }

	template<typename VECTOR>
	void mul(Grid3<T,T>&, const VECTOR&) const;

public : 
	Grid3<bool,T>		isBoundary;
	Grid3x4<T>			poisson;
};

template<typename T>
void Grid3_PoissonMatrix<T>::resize(index_type nX, index_type nY, index_type nZ, const T h, const Vector3<T>& offset) {
	isBoundary.init(nX, nY, nZ, h, h, h, offset);
	poisson.init(nX, nY, nZ);

	isBoundary.set(false);
}

template<typename T> template<typename GRID3>
void Grid3_PoissonMatrix<T>::mul(Grid3<T,T>& y, const GRID3& x) const {
	y.zero();
	compute3D_for(ApplyPoisson_Grid3<T,GRID3>(isBoundary, poisson, x, y), 1, y.nX-1, 1, y.nY-1, 1, y.nZ-1);	
}

template<typename T>
class Grid3_Preconditioner {
public : 
	typedef typename Grid3<T,T>::index_type index_type;

	Grid3_Preconditioner() {}

	virtual void buildPreconditioner(const Grid3_PoissonMatrix<T>& mat) = 0;
	virtual void applyPreconditioner(Grid3<T,T>&, const Grid3<T,T>&) const = 0;
};

template<typename T>
class Grid3_Preconditioner_Identity : public Grid3_Preconditioner<T> {
public : 
	Grid3_Preconditioner_Identity() {}

	void buildPreconditioner(const Grid3_PoissonMatrix<T>& mat) {}
	void applyPreconditioner(Grid3<T,T>& x, const Grid3<T,T>& y) const { x = y; }
};

};		// arithmetic

#endif