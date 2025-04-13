#ifndef _ARITHMETIC_GRID_PRECONDITIONER_ICCG_H_
#define _ARITHMETIC_GRID_PRECONDITIONER_ICCG_H_

#include "../../../Basis/Basis/Environments.h"
#include "Grid2_PoissonMatrix.h"

namespace arithmetic {

template<typename T>
class Grid2_Preconditioner_ICCG : public Grid2_Preconditioner<T> {
public : 
	Grid2_Preconditioner_ICCG() : matrix(NULL) {}
	Grid2_Preconditioner_ICCG(index_type nX, index_type nY) : matrix(NULL) { resize(nX,nY); }

	void resize(index_type, index_type);

	void buildPreconditioner(const Grid2_PoissonMatrix<T>& mat);
	void applyPreconditioner(Grid2<T,T>&, const Grid2<T,T>&) const;

public : 
	const Grid2_PoissonMatrix<T>* matrix;

	Grid2<T,T>				preconditioner;
	mutable Grid2<T,T>		m;
};

template<typename T> 
void Grid2_Preconditioner_ICCG<T>::resize(index_type nX, index_type nY) {
	if(preconditioner.nX == nX && preconditioner.nY == nY)
		return;

	preconditioner.init(nX, nY, 0, 0);
	m.init(nX, nY, 0, 0);
}

template<typename T>
void Grid2_Preconditioner_ICCG<T>::buildPreconditioner(const Grid2_PoissonMatrix<T>& mat) {

	index_type nX = mat.nX(), nY = mat.nY();
	resize(nX, nY);

	matrix = &mat;	

	const T micParam = 0.0f;
	real d;

	preconditioner.zero();

	for ( index_type i = 1; i < nX-1; ++i ) {
		for ( index_type j = 1; j < nY-1; ++j ) {
			if( !matrix->isBoundary(i,j) )	{

				d = matrix->poisson(i,j,0) 
					- SQ( matrix->poisson(i-1,j,1)*preconditioner(i-1,j) )
					- SQ( matrix->poisson(i,j-1,2)*preconditioner(i,j-1) )
					- micParam * ( matrix->poisson(i-1,j,1)*matrix->poisson(i-1,j,2)*SQ(preconditioner(i-1,j))
								 + matrix->poisson(i,j-1,2)*matrix->poisson(i,j-1,1)*SQ(preconditioner(i,j-1)))
								 ;

				preconditioner(i,j) = 1.0f / sqrtf((float)::abs(d)+0.001f);
				//preconditioner(i,j) = 1.0 / sqrt(::abs(d)+0.001);				
			}
		}
	}	
}

template<typename T>
void Grid2_Preconditioner_ICCG<T>::applyPreconditioner(Grid2<T,T>& y, const Grid2<T,T>& x) const {
	real d;
	index_type nX = x.nX, nY = x.nY;

	m.zero(); y.zero();

	// solve L*m = x
	for ( index_type i = 1; i < nX-1; ++i ) {
		for ( index_type j = 1; j < nY-1; ++j ) {
			if( !matrix->isBoundary(i,j) )	{
				d = x(i,j) 
					- matrix->poisson(i-1,j,1)*preconditioner(i-1,j)*m(i-1,j)
					- matrix->poisson(i,j-1,2)*preconditioner(i,j-1)*m(i,j-1);
				m(i,j) = preconditioner(i,j)*d;
			}
		}
	}
	
	// solve L'*y = m
	for ( index_type i = nX-2; i > 0; --i ) {
		for ( index_type j = nY-2; j > 0; --j ) {
			if( !matrix->isBoundary(i,j) )	{
				d = m(i,j) 
					- matrix->poisson(i,j,1)*preconditioner(i,j)*y(i+1,j)
					- matrix->poisson(i,j,2)*preconditioner(i,j)*y(i,j+1);
				y(i,j) = preconditioner(i,j)*d;
			}
		}
	}
}

};			// arithmetic

#endif