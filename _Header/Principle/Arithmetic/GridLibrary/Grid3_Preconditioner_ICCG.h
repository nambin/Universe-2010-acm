#ifndef _ARITHMETIC_GRID3_PRECONDITIONER_ICCG_H_
#define _ARITHMETIC_GRID3_PRECONDITIONER_ICCG_H_

#include "Grid3_PoissonMatrix.h"

namespace arithmetic {

template<typename T>
class Grid3_Preconditioner_ICCG : public Grid3_Preconditioner<T> {
public : 
	Grid3_Preconditioner_ICCG() : matrix(NULL) {}
	Grid3_Preconditioner_ICCG(index_type nX, index_type nY, index_type nZ) : matrix(NULL) { resize(nX,nY,nZ); }

	void resize(index_type, index_type, index_type);

	void buildPreconditioner(const Grid3_PoissonMatrix<T>& mat);
	void applyPreconditioner(Grid3<T,T>&, const Grid3<T,T>&) const;

public : 
	const Grid3_PoissonMatrix<T>* matrix;

	Grid3<T,T>				preconditioner;
	mutable Grid3<T,T>		m;
};

template<typename T> 
void Grid3_Preconditioner_ICCG<T>::resize(index_type nX, index_type nY, index_type nZ) {
	if(preconditioner.nX == nX && preconditioner.nY == nY && preconditioner.nZ == nZ)
		return;

	preconditioner.init(nX, nY, nZ, 0, 0, 0);
	m.init(nX, nY, nZ, 0, 0, 0);
}

template<typename T>
void Grid3_Preconditioner_ICCG<T>::buildPreconditioner(const Grid3_PoissonMatrix<T>& mat) {

	index_type nX = mat.nX(), nY = mat.nY(), nZ = mat.nZ();
	resize(nX, nY, nZ);

	matrix = &mat;	

	double d;

	preconditioner.zero();

	for ( index_type i = 1; i < nX-1; ++i ) {
		for ( index_type j = 1; j < nY-1; ++j ) {
			for(index_type k = 1; k < nZ-1; ++k)	{
				if( !matrix->isBoundary(i,j,k) )		{
					d = matrix->poisson(i,j,k,0) 
						- SQ( matrix->poisson(i-1,j,k,1)*preconditioner(i-1,j,k) )
						- SQ( matrix->poisson(i,j-1,k,2)*preconditioner(i,j-1,k) )
						- SQ( matrix->poisson(i,j,k-1,3)*preconditioner(i,j,k-1) );

					preconditioner(i,j,k) = 1.0f / sqrt((float)::abs(d)+0.001f);
				}
			}
		}
	}
}

template<typename T>
void Grid3_Preconditioner_ICCG<T>::applyPreconditioner(Grid3<T,T>& y, const Grid3<T,T>& x) const {

	double d;
	index_type nX = x.nX, nY = x.nY, nZ = x.nZ;

	m.zero(); y.zero();

	const bool *isBoundary_ptr = matrix->isBoundary.data;
	const double *poisson_ptr = matrix->poisson.data;
	const double *preconditioner_ptr = preconditioner.data;
	const double *x_ptr = x.data;
	double *y_ptr = y.data;
	double *m_ptr = m.data;

#define IS(i,j,k) ((k) + (nZ) * ((j) + (nY) * (i)))

	// solve L*m = x
	for ( index_type i = 1; i < nX-1; ++i ) {
		for ( index_type j = 1; j < nY-1; ++j ) {
			for( index_type k = 1; k < nZ-1; ++k )	{
				if( !isBoundary_ptr[IS(i,j,k)] )	{		
					d = x_ptr[IS(i,j,k)] 
					- poisson_ptr[4*IS(i-1,j,k)+1]*preconditioner_ptr[IS(i-1,j,k)]*m_ptr[IS(i-1,j,k)]
					- poisson_ptr[4*IS(i,j-1,k)+2]*preconditioner_ptr[IS(i,j-1,k)]*m_ptr[IS(i,j-1,k)]
					- poisson_ptr[4*IS(i,j,k-1)+3]*preconditioner_ptr[IS(i,j,k-1)]*m_ptr[IS(i,j,k-1)];
					m_ptr[IS(i,j,k)] = preconditioner_ptr[IS(i,j,k)]*d;
				}
			}
		}
	}

	// solve L'*y = m
	for ( index_type i = nX-2; i > 0; i-- ) {
		for ( index_type j = nY-2; j > 0; j-- ) {
			for( index_type k = nZ-2; k > 0; k-- )	{
				if( !isBoundary_ptr[IS(i,j,k)] )	{		
					d = m_ptr[IS(i,j,k)] 
					- poisson_ptr[4*IS(i,j,k)+1]*preconditioner_ptr[IS(i,j,k)]*y_ptr[IS(i+1,j,k)]
					- poisson_ptr[4*IS(i,j,k)+2]*preconditioner_ptr[IS(i,j,k)]*y_ptr[IS(i,j+1,k)]
					- poisson_ptr[4*IS(i,j,k)+3]*preconditioner_ptr[IS(i,j,k)]*y_ptr[IS(i,j,k+1)];
					y_ptr[IS(i,j,k)] = preconditioner_ptr[IS(i,j,k)]*d;
				}
			}
		}
	}

#undef IS
}

};			// arithmetic

#endif