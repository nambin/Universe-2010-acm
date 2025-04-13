#ifndef _VECTOR_OPERATION_FUNCTOR_H_
#define _VECTOR_OPERATION_FUNCTOR_H_

#include "../../../Basis/Basis/Functors.h"
#include "../../../Basis/Basis/Environments.h"
#include "VectorSpecialStructure.h"

namespace arithmetic {

template<typename VECX, typename R>
class VectorFunctor_SetValue {
	VECX&			x;
	const R&		val;

public : 
	typedef typename VECX::index_type index_type;

	VectorFunctor_SetValue(VECX& _x, const R& _v) : x(_x), val(_v) {}
	void operator()( const basis::blocked_range<typename VECX::index_type>& r ) const {
#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for
#else
		typename VECX::index_type my_begin = r.begin();
		typename VECX::index_type my_end = r.end();
		typename VECX::index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) 
			x[i] = val;
	}
};

template<typename VECX, typename VECA>
class VectorFunctor_Set {
	VECX&			x;
	const VECA&		a;

public : 
	typedef typename VECX::index_type index_type;

	VectorFunctor_Set(VECX& _x, const VECA& _a) : x(_x), a(_a) {}
	void operator()( const basis::blocked_range<typename VECX::index_type>& r ) const {
#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for
#else
		typename VECX::index_type my_begin = r.begin();
		typename VECX::index_type my_end = r.end();
		typename VECX::index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			x[i] = a[i];
		}
	}
};

template<typename VECX, typename VECA>
class VectorFunctor_NegativeSet {
	VECX&			x;
	const VECA&		a;

public : 
	typedef typename VECX::index_type index_type;

	VectorFunctor_NegativeSet(VECX& _x, const VECA& _a) : x(_x), a(_a) {}
	void operator()( const basis::blocked_range<typename VECX::index_type>& r ) const {
#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for
#else
		typename VECX::index_type my_begin = r.begin();
		typename VECX::index_type my_end = r.end();
		typename VECX::index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			x[i] = - a[i];
		}
	}
};

template<typename VECX, typename VECA, typename VECB>
class VectorFunctor_Add {
	VECX&		x;
	const VECA&	a;
	const VECB& b;

public :
	typedef typename VECX::index_type index_type;

	VectorFunctor_Add(VECX& _x, const VECA& _a, const VECB& _b) : x(_x), a(_a), b(_b) {}
	void operator()( const basis::blocked_range<typename VECX::index_type>& r ) const {
#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for
#else
		typename VECX::index_type my_begin = r.begin();
		typename VECX::index_type my_end = r.end();
		typename VECX::index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			x[i] = a[i] + b[i];
		}
	}
};

template<typename VECX, typename VECA, typename VECB>
class VectorFunctor_Dif {
	VECX&		x;
	const VECA&	a;
	const VECB& b;

public :
	typedef typename VECX::index_type index_type;

	VectorFunctor_Dif(VECX& _x, const VECA& _a, const VECB& _b) : x(_x), a(_a), b(_b) {}
	void operator()( const basis::blocked_range<typename VECX::index_type>& r ) const {
#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for
#else
		typename VECX::index_type my_begin = r.begin();
		typename VECX::index_type my_end = r.end();
		typename VECX::index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			x[i] = a[i] - b[i];
		}
	}
};

template<typename VECX, typename VECA>
class VectorFunctor_Icr {
	VECX&		x;
	const VECA&	a;

public :
	typedef typename VECX::index_type index_type;

	VectorFunctor_Icr(VECX& _x, const VECA& _a) : x(_x), a(_a) {}
	void operator()( const basis::blocked_range<typename VECX::index_type>& r ) const {
#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for
#else
		typename VECX::index_type my_begin = r.begin();
		typename VECX::index_type my_end = r.end();
		typename VECX::index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) 
			x[i] += a[i];
	}
};

template<typename VECX, typename VECA>
class VectorFunctor_Dcr {
	VECX&		x;
	const VECA&	a;

public :
	typedef typename VECX::index_type index_type;

	VectorFunctor_Dcr(VECX& _x, const VECA& _a) : x(_x), a(_a) {}
	void operator()( const basis::blocked_range<typename VECX::index_type>& r ) const {
#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for
#else
		typename VECX::index_type my_begin = r.begin();
		typename VECX::index_type my_end = r.end();
		typename VECX::index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) 
			x[i] -= a[i];
	}
};

template<typename VECX, typename VECA, typename R>
class VectorFunctor_Mul {
	VECX&		x;
	const VECA&	a;
	const R		c;

public :
	typedef typename VECX::index_type index_type;

	VectorFunctor_Mul(VECX& _x, const VECA& _a, const R _c) : x(_x), a(_a), c(_c) {}
	void operator()( const basis::blocked_range<typename VECX::index_type>& r ) const {
#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for
#else
		typename VECX::index_type my_begin = r.begin();
		typename VECX::index_type my_end = r.end();
		typename VECX::index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) 
			x[i] = (a[i] * c);
	}
};

template<typename VECX, typename R>
class VectorFunctor_Scale {
	VECX&		x;
	const R		c;

public :
	typedef typename VECX::index_type index_type;

	VectorFunctor_Scale(VECX& _x, const R _c) : x(_x), c(_c) {}
	void operator()( const basis::blocked_range<typename VECX::index_type>& r ) const {
#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for
#else
		typename VECX::index_type my_begin = r.begin();
		typename VECX::index_type my_end = r.end();
		typename VECX::index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) 
			x[i] *= c;
	}
};

template<typename VECX, typename VECY>
class VectorFunctor_Dot {
	typedef typename VECX::value_type T;

	const VECX&	x;
	const VECY&	y;
	mutable T	val;

public :
	typedef typename VECX::index_type index_type;

	VectorFunctor_Dot(const VECX& _x, const VECY& _y) : x(_x), y(_y), val(0) {}	
	const T getValue() const { return val; }

	void operator()( const basis::blocked_range<typename VECX::index_type>& r ) const {
		T _val = T(0.0);

#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for reduction(+:_val) 
#else
		typename VECX::index_type my_begin = r.begin();
		typename VECX::index_type my_end = r.end();
		typename VECX::index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			_val += (x[i] * y[i]);
		}
		val = _val;
	}	
};

template<typename VECX>
class VectorFunctor_Dot_Directly {
	typedef typename VECX::value_type T;

	const T* const	x;
	const T* const	y;
	mutable T		val;

public :
	typedef typename VECX::index_type index_type;

	VectorFunctor_Dot_Directly(const VECX& _x, const VECX& _y) : x(_x.getptr()), y(_y.getptr()), val(0) {}	
	const T getValue() const { return val; }

	void operator()( const basis::blocked_range<typename VECX::index_type>& r ) const {
		T _val = T(0.0);

#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for reduction(+:_val) 
#else
		typename VECX::index_type my_begin = r.begin();
		typename VECX::index_type my_end = r.end();
		typename VECX::index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			_val += (x[i] * y[i]);
		}
		val = _val;
	}	
};

template<typename VECX, typename VECA, typename R>
class VectorFunctor_IcrScaled {
	VECX&		x;
	const VECA&	a;
	const R		c;

public :
	typedef typename VECX::index_type index_type;

	VectorFunctor_IcrScaled(VECX& _x, const VECA& _a, const R _c) : x(_x), a(_a), c(_c) {}
	void operator()( const basis::blocked_range<typename VECX::index_type>& r ) const {
#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for
#else
		typename VECX::index_type my_begin = r.begin();
		typename VECX::index_type my_end = r.end();
		typename VECX::index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) 
			x[i] += (c * a[i]);
	}
};

template<typename VECX>
class VectorFunctor_IcrScaled_Directly {
	typedef typename VECX::value_type R;

	R*				x;
	const R* const	a;
	const R			c;

public :
	typedef typename VECX::index_type index_type;

	VectorFunctor_IcrScaled_Directly(VECX& _x, const VECX& _a, const R _c) : x(_x.getptr()), a(_a.getptr()), c(_c) {}
	void operator()( const basis::blocked_range<typename VECX::index_type>& r ) const {
#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for
#else
		typename VECX::index_type my_begin = r.begin();
		typename VECX::index_type my_end = r.end();
		typename VECX::index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) 
			x[i] += (c * a[i]);
	}
};

template<typename VECX, typename VECA, typename R>
class VectorFunctor_ScaledIcr {
	VECX&		x;
	const VECA&	a;
	const R		c;

public :
	typedef typename VECX::index_type index_type;

	VectorFunctor_ScaledIcr(VECX& _x, const VECA& _a, const R _c) : x(_x), a(_a), c(_c) {}
	void operator()( const basis::blocked_range<typename VECX::index_type>& r ) const {
#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for
#else
		typename VECX::index_type my_begin = r.begin();
		typename VECX::index_type my_end = r.end();
		typename VECX::index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) 
			x[i] = a[i] + c * x[i];
	}
};

template<typename VECX>
class VectorFunctor_ScaledIcr_Directly {
	typedef typename VECX::value_type R;

	R*				x;
	const R* const	a;
	const R			c;

public :
	typedef typename VECX::index_type index_type;

	VectorFunctor_ScaledIcr_Directly(VECX& _x, const VECX& _a, const R _c) : x(_x.getptr()), a(_a.getptr()), c(_c) {}
	void operator()( const basis::blocked_range<typename VECX::index_type>& r ) const {
#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for
#else
		typename VECX::index_type my_begin = r.begin();
		typename VECX::index_type my_end = r.end();
		typename VECX::index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) 
			x[i] = a[i] + c * x[i];
	}
};

template<typename VECX>
class VectorFunctor_Sum {
	typedef typename VECX::value_type T;

	const VECX&	x;
	mutable T	val;

public :
	typedef typename VECX::index_type index_type;

	VectorFunctor_Sum(const VECX& _x) : x(_x), val(0) {}	
	const T getValue() const { return val; }

	void operator()( const basis::blocked_range<typename VECX::index_type>& r ) const {
		T _val = T(0.0);

#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for reduction(+:_val) 
#else
		typename VECX::index_type my_begin = r.begin();
		typename VECX::index_type my_end = r.end();
		typename VECX::index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			_val += x[i];
		}
		val = _val;
	}	
};

template<typename VECX, typename VECY>
class VectorFunctor_SqrLength {
	typedef typename VECX::value_type T;

	const VECX&	x;
	const VECY&	y;
	mutable T	val;

public :
	typedef typename VECX::index_type index_type;

	VectorFunctor_SqrLength(const VECX& _x, const VECY& _y) : x(_x), y(_y), val(0) {}
	const T getValue() const { return val; }

	void operator()( const basis::blocked_range<typename VECX::index_type>& r ) const {
		T _val = T(0.0);

#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for reduction(+:_val) 
#else
		typename VECX::index_type my_begin = r.begin();
		typename VECX::index_type my_end = r.end();
		typename VECX::index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			_val += (x[i] - y[i]) * (x[i] - y[i]);	
		}
		val = _val;
	}	
};

template<typename VECX>
class VectorFunctor_L1Norm {
	typedef typename VECX::value_type T;

	const T* const	x;
	mutable T		val;

public :
	typedef typename VECX::index_type index_type;

	VectorFunctor_L1Norm(const VECX& _x) : x(_x.getptr()), val(0) {}	
	const T getValue() const { return val; }

	void operator()( const basis::blocked_range<typename VECX::index_type>& r ) const {
		T _val = T(0.0);

#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for reduction(+:_val) 
#else
		typename VECX::index_type my_begin = r.begin();
		typename VECX::index_type my_end = r.end();
		typename VECX::index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			_val += ::abs(x[i]);
		}
		val = _val;
	}	
};

template<typename VECX>
class VectorFunctor_L2Norm {
	typedef typename VECX::value_type T;

	const T* const	x;
	mutable T		val;

public :
	typedef typename VECX::index_type index_type;

	VectorFunctor_L2Norm(const VECX& _x) : x(_x.getptr()), val(0) {}	
	const T getValue() const { return sqrt(val); }

	void operator()( const basis::blocked_range<typename VECX::index_type>& r ) const {
		T _val = T(0.0);

#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for reduction(+:_val) 
#else
		typename VECX::index_type my_begin = r.begin();
		typename VECX::index_type my_end = r.end();
		typename VECX::index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			_val += (x[i] * x[i]);
		}
		val = _val;
	}	
};

template<typename VECX>
class VectorFunctor_InfNorm {
	typedef typename VECX::value_type T;

	const T* const	x;
	mutable T		val;

public :
	typedef typename VECX::index_type index_type;

	VectorFunctor_InfNorm(const VECX& _x) : x(_x.getptr()), val(0) {}	
	const T getValue() const { return val; }

	void operator()( const basis::blocked_range<typename VECX::index_type>& r ) const {

		val = 0;

#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for
		for( i = my_begin; i < my_end; ++i ) {
			if(val < ::abs(x[i])) {
#pragma omp critical
				{	
					// ret could have changed by another thread
					// so, compare again
					val = (::abs(x[i]) > val) ? ::abs(x[i]) : val;
				}
			}
		}
#else
		for( typename VECX::index_type i = r.begin(); i < r.end(); ++i ) {
			if(val < ::abs(x[i])) val = ::abs(x[i]);
		}
#endif
	}	
};

template<typename VECX, typename VECA, typename VECB>
class VectorFunctor_MaxVec {
	VECX&		x;
	const VECA&	a;
	const VECB&	b;

public :
	typedef typename VECX::index_type index_type;

	VectorFunctor_MaxVec(VECX& _x, const VECA& _a, const VECB& _b) : x(_x), a(_a), b(_b) {}
	void operator()( const basis::blocked_range<typename VECX::index_type>& r ) const {
#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for
#else
		typename VECX::index_type my_begin = r.begin();
		typename VECX::index_type my_end = r.end();
		typename VECX::index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) 
			x[i] = MAX(a[i], b[i]);
	}
};

template<typename VECX, typename VECA, typename VECB>
class VectorFunctor_MinVec {
	VECX&		x;
	const VECA&	a;
	const VECB&	b;

public :
	typedef typename VECX::index_type index_type;

	VectorFunctor_MinVec(VECX& _x, const VECA& _a, const VECB& _b) : x(_x), a(_a), b(_b) {}
	void operator()( const basis::blocked_range<typename VECX::index_type>& r ) const {
#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for
#else
		typename VECX::index_type my_begin = r.begin();
		typename VECX::index_type my_end = r.end();
		typename VECX::index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) 
			x[i] = MIN(a[i], b[i]);
	}
};

};		// arithmetic

#endif