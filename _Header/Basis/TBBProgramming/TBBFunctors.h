#pragma once

#include <cfloat>
#include <cmath>
#include <tbb/blocked_range.h>
#include <tbb/blocked_range2d.h>
#include <tbb/blocked_range3d.h>

namespace basis {

template<typename FUNCTOR>
class TBBFunctor_MinIndex {
	typedef typename FUNCTOR::ret_type ret_type;
	typedef typename FUNCTOR::index_type index_type;

	const FUNCTOR& func;

public:
	mutable ret_type		value;
	mutable index_type		index;

	void operator()( const tbb::blocked_range<index_type>& r ) const {
		for( index_type i = r.begin(); i != r.end( ); ++i ) {
			ret_type v = func(i);
			if( v < value ) {
				value = v;
				index = i;
			}
		}
	}

	TBBFunctor_MinIndex( const FUNCTOR& _func ) : value(ret_type(FLT_MAX)), index(0), func(_func) {}
	TBBFunctor_MinIndex( TBBFunctor_MinIndex& x, tbb::split ) : value(ret_type(FLT_MAX)), index(0), func(x.func) {}
	void join( const TBBFunctor_MinIndex& y ) {
		if( y.value < value ) {
			value = y.value;
			index = y.index;
		}
	}	
};

template<typename FUNCTOR>
class TBBFunctor_MinIndex2d {
	typedef typename FUNCTOR::ret_type ret_type;
	typedef typename FUNCTOR::index_type index_type;

	const FUNCTOR& func;

public:
	mutable ret_type		value;
	mutable index_type		index0;
	mutable index_type		index1;

	void operator()( const tbb::blocked_range2d<index_type>& r ) const {
		for( index_type i = r.rows().begin(); i != r.rows().end( ); ++i ) {
			for( index_type j = r.cols().begin(); j != r.cols().end( ); ++j ) {
				ret_type v = func(i,j);
				if( v < value ) {
					value = v;
					index0 = i;
					index1 = j;
				}
			}			
		}
	}

	TBBFunctor_MinIndex2d( const FUNCTOR& _func ) : value(ret_type(FLT_MAX)), index0(0), index1(0), func(_func) {}
	TBBFunctor_MinIndex2d( TBBFunctor_MinIndex2d& x, tbb::split ) : value(ret_type(FLT_MAX)), index0(0), index1(0), func(x.func) {}
	void join( const TBBFunctor_MinIndex2d& y ) {
		if( y.value < value ) {
			value = y.value;
			index0 = y.index0;
			index1 = y.index1;
		}
	}	
};

template<typename FUNCTOR>
class TBBFunctor_MinIndex3d {
	typedef typename FUNCTOR::ret_type ret_type;
	typedef typename FUNCTOR::index_type index_type;

	const FUNCTOR& func;

public:
	mutable ret_type		value;
	mutable index_type		index0;
	mutable index_type		index1;
	mutable index_type		index2;

	void operator()( const tbb::blocked_range3d<index_type>& r ) const {
		for( index_type i = r.pages().begin(); i != r.pages().end( ); ++i ) {
			for( index_type j = r.rows().begin(); j != r.rows().end( ); ++j ) {
				for( index_type k = r.cols().begin(); k != r.cols().end( ); ++k ) {
					ret_type v = func(i,j,k);
					if( v < value ) {
						value = v;
						index0 = i;
						index1 = j;
						index2 = k;
					}
				}			
			}
		}
	}

	TBBFunctor_MinIndex3d( const FUNCTOR& _func ) : value(ret_type(FLT_MAX)), index0(0), index1(0), index2(0), func(_func) {}
	TBBFunctor_MinIndex3d( TBBFunctor_MinIndex3d& x, tbb::split ) : value(ret_type(FLT_MAX)), index0(0), index1(0), index2(0), func(x.func) {}
	void join( const TBBFunctor_MinIndex3d& y ) {
		if( y.value < value ) {
			value = y.value;
			index0 = y.index0;
			index1 = y.index1;
			index2 = y.index2;
		}
	}	
};

};		// basis