#ifndef _ARITHMETIC_GRID3_POISSON_FUNCTOR_H_
#define _ARITHMETIC_GRID3_POISSON_FUNCTOR_H_

#include "../../../Basis/Basis/Environments.h"
#include "GridLibrary.h"

namespace arithmetic {

template<typename T, typename VECTOR>
class ApplyPoisson_Grid3 {
public : 
	typedef typename Grid3<T,T>::index_type index_type;

	const Grid3<bool,T>&	isBoundary;
	const Grid3x4<T>&		poisson;
	const VECTOR&			x;
	Grid3<T,T>&				y;

	ApplyPoisson_Grid3(const Grid3<bool,T>& _isBdry, const Grid3x4<T>& _poi, const VECTOR& _x, Grid3<T,T>& _y) 
	: isBoundary(_isBdry), poisson(_poi), x(_x), y(_y) {}

	void operator()( const basis::blocked_range3d<index_type>& r ) const {
		index_type nY = y.nY;
		index_type nZ = y.nZ;

#ifdef USE_OPENMP
		int my_begin = int(r.pages().begin());
		int my_end = int(r.pages().end());
		int i;

#pragma omp parallel for
#else
		index_type my_begin = r.pages().begin();
		index_type my_end = r.pages().end();
		index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			for( index_type j=r.rows().begin(); j!=r.rows().end( ); ++j ) {
				for( index_type k=r.cols().begin(); k!=r.cols().end( ); ++k ) {	

#ifdef USE_GRID_INLINE_FUNC
					if( !isBoundary(i,j,k) )	{
						y(i,j,k)	= poisson(i,j,k,0)   * x(i,j,k)
									+ poisson(i,j,k-1,3) * x(i,j,k-1)
									+ poisson(i,j,k,3)   * x(i,j,k+1)
									+ poisson(i-1,j,k,1) * x(i-1,j,k)
									+ poisson(i,j,k,1)   * x(i+1,j,k)
									+ poisson(i,j-1,k,2) * x(i,j-1,k)
									+ poisson(i,j,k,2)   * x(i,j+1,k);
					}

#else		// #ifdef USE_GRID_INLINE_FUNC

#define IS(i,j,k) ((k) + (nZ) * ((j) + (nY) * (i)))
					if( !isBoundary[IS(i,j,k)] )	{
						y[IS(i,j,k)]	= poisson[4*IS(i,j,k)]     * x[IS(i,j,k)]
										+ poisson[4*IS(i,j,k-1)+3] * x[IS(i,j,k-1)]
										+ poisson[4*IS(i,j,k)+3]   * x[IS(i,j,k+1)]
										+ poisson[4*IS(i-1,j,k)+1] * x[IS(i-1,j,k)]
										+ poisson[4*IS(i,j,k)+1]   * x[IS(i+1,j,k)]
										+ poisson[4*IS(i,j-1,k)+2] * x[IS(i,j-1,k)]
										+ poisson[4*IS(i,j,k)+2]   * x[IS(i,j+1,k)];
#undef IS

#endif		// #ifdef USE_GRID_INLINE_FUNC
					}					
				}		
			}
		}


	}
};

};		// arithmetic

#endif