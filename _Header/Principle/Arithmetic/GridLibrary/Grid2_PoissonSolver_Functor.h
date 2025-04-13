#ifndef _ARITHMETIC_GRID2_POISSON_FUNCTOR_H_
#define _ARITHMETIC_GRID2_POISSON_FUNCTOR_H_

#include "GridLibrary.h"

namespace arithmetic {

template<typename R, typename VECTOR>
class ApplyPoisson_Grid2 {
public : 
	typedef typename Grid2<R,R>::index_type index_type;

	const Grid2<bool,R>&	isBoundary;
	const Grid2x3<R>&		poisson;
	const VECTOR&			x;
	Grid2<R,R>&				y;

	ApplyPoisson_Grid2(const Grid2<bool,R>& _isBdry, const Grid2x3<R>& _poi, const VECTOR& _x, Grid2<R,R>& _y) 
	: isBoundary(_isBdry), poisson(_poi), x(_x), y(_y) {}

	void operator()( const basis::blocked_range2d<index_type>& r ) const {


	}
};

};		// arithmetic

#endif