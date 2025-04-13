#ifndef _ARITHMETIC_VECTOR_SPECIAL_STRUCTURE_H_
#define _ARITHMETIC_VECTOR_SPECIAL_STRUCTURE_H_

#include "../../../Basis/Basis/Environments.h"

namespace arithmetic {

template<typename VECTOR, typename REAL_TYPE>
class Scaled {
public : 
	typedef typename VECTOR::value_type value_type;
	typedef typename VECTOR::index_type	index_type;

	Scaled(const VECTOR& _v, REAL_TYPE _s) : v(_v), s(_s) {}

	index_type size() const { return v.size(); }

	const value_type operator[](index_type i) const { return s * v[i]; }	

	const value_type operator()(index_type i) const { return s * v(i); }
	const value_type operator()(index_type i, index_type j) const { return s * v(i,j); }
	const value_type operator()(index_type i, index_type j, index_type k) const { return s * v(i,j,k); }

protected : 
	const VECTOR&		v;
	REAL_TYPE			s;
};

template<typename VECTOR, typename REAL_TYPE>
Scaled<VECTOR, REAL_TYPE> scaled(const VECTOR& vec, REAL_TYPE c) {
	return Scaled<VECTOR, REAL_TYPE>(vec, c);
}

};		// arithmetic

#endif