#pragma once

#include "TBBFunctors.h"
#include <tbb/parallel_reduce.h>

namespace basis {

template<typename FUNCTOR>
typename FUNCTOR::index_type TBBFunction_MinIndex( const FUNCTOR &func, typename FUNCTOR::index_type s, typename FUNCTOR::index_type n, 
												   typename tbb::blocked_range<typename FUNCTOR::index_type>::size_type grainSize, bool useAutoPartitioner = true ) 
{
	TBBFunctor_MinIndex<FUNCTOR> minIdx(func);
	if(useAutoPartitioner)	tbb::parallel_reduce(blocked_range<typename FUNCTOR::index_type>(s,n), minIdx, auto_partitioner() );
	else					tbb::parallel_reduce(blocked_range<typename FUNCTOR::index_type>(s,n,grainSize), minIdx );
	return minIdx.index;
}

};		// basis