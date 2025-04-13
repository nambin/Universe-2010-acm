#ifndef _ARITHMETIC_ARRAY_FUNCTOR_H_
#define _ARITHMETIC_ARRAY_FUNCTOR_H_

#include "../../../Basis/Basis/Environments.h"

namespace arithmetic {

template<typename T> class ArrayNegativeAbs {
	const T* const ptr;
public : 
	typedef T ret_type;

	ArrayNegativeAbs(const T* const _ptr) : ptr(_ptr) {}
	inline T operator()(size_t i) const { return - ::abs(ptr[i]); }
};

template<typename T> class ArrayFunctor {
	const T* const ptr;
public : 
	ArrayFunctor(const T* const _ptr) : ptr(_ptr) {}
	inline T operator()(size_t i) const { return ptr[i]; }
};

template<typename T> class ArraySetFunctor {
	T* ptr;
	T val;

public : 
	ArraySetFunctor(T *_ptr, T _val) : ptr(_ptr), val(_val) {}
	template<typename INDEX> void operator()( const basis::blocked_range<INDEX>& r ) const {
#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for
#else
		INDEX my_begin = r.begin();
		INDEX my_end = r.end();
		INDEX i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			ptr[i] = val;
		}
	}
};

template<typename T> class ArrayDotFunctor {
	const T* const ptr0;
	const T* const ptr1;
	mutable T val;

public : 	
	template<typename INDEX> void operator()( const basis::blocked_range<INDEX>& r ) const {
		T _val = T(0.0);

#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for reduction(+:_val) 
#else
		INDEX my_begin = r.begin();
		INDEX my_end = r.end();
		INDEX i;
#endif
		for( i = my_begin; i < my_end; ++i ) {		
			_val += (ptr0[i] * ptr1[i]);
		}
		val = _val;
	}

	ArrayDotFunctor(const T* const _ptr0, const T* const _ptr1) : val(0), ptr0(_ptr0), ptr1(_ptr1) {}	
	const T getValue() const { return val; }
};

template<typename T> class ArrayIncrement {
	T* ptr;
	T scale;
	const T* const ptr1;	

public : 
	ArrayIncrement(T *_ptr, T _scale, T *_ptr1) : ptr(_ptr), scale(_scale), ptr1(_ptr1) {}
	template<typename INDEX> void operator()( const basis::blocked_range<INDEX>& r ) const {
#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for
#else
		INDEX my_begin = r.begin();
		INDEX my_end = r.end();
		INDEX i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			ptr[i] += scale * ptr1[i];
		}
	}	
};

template<typename T> class ArrayScaleAndIncrement {
	T* ptr;
	T scale;
	const T* const ptr1;	

public : 
	ArrayScaleAndIncrement(T *_ptr, T _scale, T *_ptr1) : ptr(_ptr), scale(_scale), ptr1(_ptr1) {}
	template<typename INDEX> void operator()( const basis::blocked_range<INDEX>& r ) const {
#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for
#else
		INDEX my_begin = r.begin();
		INDEX my_end = r.end();
		INDEX i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			ptr[i] = ptr1[i] + scale * ptr[i];			
		}
	}	
};

};		// arithmetic

#endif