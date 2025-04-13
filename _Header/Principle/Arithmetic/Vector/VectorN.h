#ifndef _ARITHMETIC_VECTORN_H_
#define _ARITHMETIC_VECTORN_H_

#include "../../../Basis/Basis/BasicLibrary.h"
#include "../MathLibrary/MatrixVectorOperation.h"
#include "VectorOperation.h"

namespace arithmetic {

#ifdef _DEBUG
#define __ASSERT(c) assert(c)
#else
#define __ASSERT(c) {}
#endif

//////////////////////////////////////////////////////////////////////////
//	VectorN class definition
//////////////////////////////////////////////////////////////////////////
template <typename T> class VectorN {
public:
	typedef std::size_t			index_type;
	typedef index_type*			index_pointer ;
	typedef const index_type*	index_const_pointer ;
	typedef index_type&			index_reference ;
	typedef const index_type&	index_const_reference ;

	typedef T					value_type;
	typedef value_type*			value_pointer ;       
	typedef const value_type*	value_const_pointer ; 
	typedef value_type&			value_reference ;     
	typedef const value_type&	value_const_reference;

	//////////////////////////////////////////////////////////////////////////
	// Constructors
	explicit VectorN(index_type N = 0, T d = 0);

	VectorN(const VectorN &src);	
	VectorN(index_type N, const T *src);	

	~VectorN();

	void release();

	//////////////////////////////////////////////////////////////////////////
	// Index operators
	inline const T  get(index_type i) const;		
	
	inline T&		operator[](index_type i);		
	inline const T	operator[](index_type i) const;	

	inline T&		operator()(index_type i);		
	inline const T  operator()(index_type i) const;	

	//////////////////////////////////////////////////////////////////////////
	// Size
	inline const index_type dim()	const { return _dim; }
	inline const index_type size()	const { return _dim; }	

	//////////////////////////////////////////////////////////////////////////
	// Resize 
	void resize( index_type n );

	//////////////////////////////////////////////////////////////////////////
	// Assignment and set
	void zero();

	void set(T d);			
	inline void set(index_type i, T d);	

	void set(const VectorN &src);
	void set(const T *src);		

	const VectorN &operator=(T d);				
	const VectorN &operator=(const VectorN &src);
	const VectorN &operator=(const T *src);	
	const VectorN &operator=(const Vector_neg_V<VectorN<T>> & op);
	const VectorN &operator=(const Vector_V_add_V<VectorN<T>> & op);
	const VectorN &operator=(const Vector_V_dif_V<VectorN<T>> & op);
	const VectorN &operator=(const Vector_V_mul_S<VectorN<T>,T> & op);

	template<typename MATRIX> const VectorN &operator=(const Vector_M_mul_V<MATRIX, VectorN<T>> &op);
	template<typename MATRIX> const VectorN &operator=(const Vector_S_mul_M_mul_V<T, MATRIX, VectorN<T>> &op);
	template<typename MATRIX> const VectorN &operator=(const Vector_V_sum_M_mul_V<MATRIX, VectorN<T>> &op);
	template<typename MATRIX> const VectorN &operator=(const Vector_V_sub_M_mul_V<MATRIX, VectorN<T>> &op);
	template<typename MATRIX> const VectorN &operator=(const Vector_V_sum_S_mul_M_mul_V<T, MATRIX, VectorN<T>> &op);
	template<typename MATRIX> const VectorN &operator=(const Vector_V_sub_S_mul_M_mul_V<T, MATRIX, VectorN<T>> &op);

	//////////////////////////////////////////////////////////////////////////
	// Manipulation
	const VectorN &normalize();
	const T mag() const;
	const T sqrmag() const;

	//////////////////////////////////////////////////////////////////////////
	// Get Pointer
	inline const T* const getptr() const { return data; }

	//////////////////////////////////////////////////////////////////////////
	// Comparison operators
	const bool operator==(const VectorN &) const;	
	const bool operator!=(const VectorN &) const;	

	const bool operator==(const T d) const;	
	const bool operator!=(const T d) const;	

	//////////////////////////////////////////////////////////////////////////
	// In place arithmetic
	const VectorN &operator+=(T d);	
	const VectorN &operator-=(T d);	
	const VectorN &operator*=(T d);	
	const VectorN &operator/=(T d);	

	//////////////////////////////////////////////////////////////////////////
	// In place arithmetic; componentwise operations
	const VectorN &operator+=(const VectorN &src);		
	const VectorN &operator-=(const VectorN &src);		
	const VectorN &operator*=(const VectorN &src);		
	const VectorN &operator/=(const VectorN &src);		

	const VectorN &operator+=(const Vector_V_mul_S<VectorN<T>,T> & op);		
	const VectorN &operator-=(const Vector_V_mul_S<VectorN<T>,T> & op);		
	const VectorN &operator*=(const Vector_V_mul_S<VectorN<T>,T> & op);	
	const VectorN &operator/=(const Vector_V_mul_S<VectorN<T>,T> & op);	

	//////////////////////////////////////////////////////////////////////////
	// IO
	template <typename T> friend std::ostream& operator<<(std::ostream &os, const VectorN<T> &v);	

	//////////////////////////////////////////////////////////////////////////
	// Friends Operators
	template <typename T> friend void neg(VectorN<T> &ret, const VectorN<T> &a);
	template <typename T> friend void add(VectorN<T> &ret, const VectorN<T> &a, const VectorN<T> &b);
	template <typename T> friend void dif(VectorN<T> &ret, const VectorN<T> &a, const VectorN<T> &b);

	template <typename T> friend void mul(VectorN<T> &ret, const VectorN<T> &v, const T c);
	template <typename T> friend void mul(VectorN<T> &ret, const T c, const VectorN<T> &v);

	template <typename T> friend void div(VectorN<T> &ret, const VectorN<T> &v, const T c);
		
	template <typename T> friend void unitVector(VectorN<T> &ret, const VectorN<T> &a);
	template <typename T> friend const T normalize(VectorN<T> &a);

	template <typename T> friend const T l1Norm(const VectorN<T> &a);
	template <typename T> friend const T l2Norm(const VectorN<T> &a);
	template <typename T> friend const T lInfNorm(const VectorN<T> &a);
	template <typename T> friend const T mag(const VectorN<T> &a);
	template <typename T> friend const T sqrmag(const VectorN<T> &a);
	template <typename T> friend const T sum(const VectorN<T> &a);	

	template <typename T> friend const T dot(const VectorN<T> &a,const VectorN<T> &b);
	template <typename T> friend const T length(const VectorN<T> &a, const VectorN<T> &b);
	template <typename T> friend const T sqrlength(const VectorN<T> &a, const VectorN<T> &b);

	template <typename T> friend void maxVec(VectorN<T> &ret,const VectorN<T> &a,const VectorN<T> &b);
	template <typename T> friend void minVec(VectorN<T> &ret,const VectorN<T> &a,const VectorN<T> &b);

	template <typename T> friend const Vector_neg_V<VectorN<T>> operator-(const VectorN<T> &a);
	template <typename T> friend const Vector_V_add_V<VectorN<T>> operator+(const VectorN<T> &a, const VectorN<T> &b);
	template <typename T> friend const Vector_V_dif_V<VectorN<T>> operator-(const VectorN<T> &a, const VectorN<T> &b);
	template <typename T> friend const Vector_V_mul_S<VectorN<T>,T> operator*(const VectorN<T> &a, const T &c);
	template <typename T> friend const Vector_V_mul_S<VectorN<T>,T> operator*(const T &c, const VectorN<T> &a);
	template <typename T> friend const Vector_V_mul_S<VectorN<T>,T> operator/(const VectorN<T> &a, const T &c);

private : 
	void init(index_type N);
	
private :
	T*			data;
	index_type	_dim;
};

//////////////////////////////////////////////////////////////////////////
// Constructors
template <typename T> VectorN<T>::VectorN( index_type N, T d ) : data(NULL)			 { init(N); set(d);				}
template <typename T> VectorN<T>::VectorN( const VectorN &src ) : data(NULL)		 { init(src.dim()); set(src);	}
template <typename T> VectorN<T>::VectorN( index_type N, const T *src ) : data(NULL) { init(N);	set(src);			}

template <typename T> VectorN<T>::~VectorN() { release(); }

template <typename T> void VectorN<T>::release() {
	if(data)
		delete [] data;
	data = NULL;
	_dim = 0;
}

template <typename T> void VectorN<T>::init(index_type N) {
	__ASSERT(N >= 0);

	_dim = N;

	if(N == 0)	data = NULL;		
	else		data = new T[_dim];	
}

//////////////////////////////////////////////////////////////////////////
// Index operators
template <typename T> inline const T VectorN<T>::get(index_type i) const		{ __ASSERT( 0 <= i && i < _dim ); return data[i]; }

template <typename T> inline T &VectorN<T>::operator[](index_type i)			{ __ASSERT( 0 <= i && i < _dim ); return data[i]; }
template <typename T> inline T &VectorN<T>::operator()(index_type i)			{ __ASSERT( 0 <= i && i < _dim ); return data[i]; }
template <typename T> inline const T VectorN<T>::operator[](index_type i) const { __ASSERT( 0 <= i && i < _dim ); return data[i]; }
template <typename T> inline const T VectorN<T>::operator()(index_type i) const { __ASSERT( 0 <= i && i < _dim ); return data[i]; }

//////////////////////////////////////////////////////////////////////////
// Resize
template <typename T> void VectorN<T>::resize( index_type n ) {
	__ASSERT(n > 0);

	release();

	init(n);
	zero();
}

//////////////////////////////////////////////////////////////////////////
// Assignment and set
template <typename T> void VectorN<T>::zero()	{ return set(T(0));								}
template <typename T> void VectorN<T>::set(T d) { for(index_type i=0;i<_dim;++i) data[i] = d;	}

template <typename T> inline void VectorN<T>::set(index_type i, T d) { __ASSERT(0<=i && i<_dim); data[i] = d; }
template <typename T> void VectorN<T>::set(const VectorN &src)		 { __ASSERT(_dim==src.dim()); for(index_type i=0;i<_dim;++i) data[i] = src[i]; }
template <typename T> void VectorN<T>::set(const T *src)			 {							  for(index_type i=0;i<_dim;++i) data[i] = src[i]; }

template <typename T> const VectorN<T> &VectorN<T>::operator=(T d)				  { set(d);   return (*this); }
template <typename T> const VectorN<T> &VectorN<T>::operator=(const VectorN &src) {	set(src); return (*this); }
template <typename T> const VectorN<T> &VectorN<T>::operator=(const T *src)		  { set(src); return (*this); }

//////////////////////////////////////////////////////////////////////////
// Vector Operation
template <typename T> const VectorN<T> &VectorN<T>::operator=(const Vector_neg_V<VectorN<T>> & op)		{ arithmetic::neg(*this, op.V);			return (*this); }
template <typename T> const VectorN<T> &VectorN<T>::operator=(const Vector_V_add_V<VectorN<T>> & op)	{ arithmetic::add(*this, op.V1, op.V2); return (*this); }
template <typename T> const VectorN<T> &VectorN<T>::operator=(const Vector_V_dif_V<VectorN<T>> & op)	{ arithmetic::dif(*this, op.V1, op.V2); return (*this); }
template <typename T> const VectorN<T> &VectorN<T>::operator=(const Vector_V_mul_S<VectorN<T>,T> & op)	{ arithmetic::mul(*this, op.V, op.s);	return (*this); }

//////////////////////////////////////////////////////////////////////////
// Matrix/Vector Operation
template<typename T> template<typename MATRIX> const VectorN<T> &VectorN<T>::operator=(const Vector_M_mul_V<MATRIX, VectorN<T>> &op) {
	mul(*this, op.M, op.a);
	return (*this);
}

template<typename T> template<typename MATRIX> const VectorN<T> &VectorN<T>::operator=(const Vector_S_mul_M_mul_V<T, MATRIX, VectorN<T>> &op) {
	mul(*this, op.M, op.a);
	(*this) *= op.s;
	return (*this);
}

template<typename T> template<typename MATRIX> const VectorN<T> &VectorN<T>::operator=(const Vector_V_sum_M_mul_V<MATRIX, VectorN<T>> &op) {
	if(&op.a != this) {
		mul(*this, op.M, op.b);
		arithmetic::add(*this, op.a, *this);
	}
	else {
		VectorN<T> _tmp(op.a);
		mul(*this, op.M, op.b);
		arithmetic::add(*this, _tmp, *this);
	}
	
	return (*this);
}

template<typename T> template<typename MATRIX> const VectorN<T> &VectorN<T>::operator=(const Vector_V_sub_M_mul_V<MATRIX, VectorN<T>> &op) {
	if(&op.a != this) {
		mul(*this, op.M, op.b);
		arithmetic::dif(*this, op.a, *this);
	}
	else {
		VectorN<T> _tmp(op.a);
		mul(*this, op.M, op.b);
		arithmetic::dif(*this, _tmp, *this);
	}
	
	return (*this);
}

template<typename T> template<typename MATRIX> const VectorN<T> &VectorN<T>::operator=(const Vector_V_sum_S_mul_M_mul_V<T, MATRIX, VectorN<T>> &op) {
	if(&op.a != this) {
		mul(*this, op.M, op.b);
		(*this) *= op.s;
		arithmetic::add(*this, op.a, *this);
	}
	else {
		VectorN<T> _tmp(op.a);
		mul(*this, op.M, op.b);
		(*this) *= op.s;
		arithmetic::add(*this, _tmp, *this);
	}
	
	return (*this);
}

template<typename T> template<typename MATRIX> const VectorN<T> &VectorN<T>::operator=(const Vector_V_sub_S_mul_M_mul_V<T, MATRIX, VectorN<T>> &op) {
	if(&op.a != this) {
		mul(*this, op.M, op.b);
		(*this) *= op.s;
		arithmetic::dif(*this, op.a, *this);
	}
	else {
		VectorN<T> _tmp(op.a);
		mul(*this, op.M, op.b);
		(*this) *= op.s;
		arithmetic::dif(*this, _tmp, *this);
	}

	return (*this);
}

//////////////////////////////////////////////////////////////////////////
// Manipulation
template <typename T> const VectorN<T> &VectorN<T>::normalize() { arithmetic::normalize(*this);	return (*this); }
template <typename T> const T VectorN<T>::mag() const			{ return arithmetic::mag(*this);				}
template <typename T> const T VectorN<T>::sqrmag() const		{ return arithmetic::sqrmag(*this);				}

//////////////////////////////////////////////////////////////////////////
// comparison Operator
template <typename T> const bool VectorN<T>::operator==( const VectorN &src ) const {
	__ASSERT( _dim == src.dim() );
	for(index_type i = 0; i < _dim; ++i) {
		if(!IsAlmostSame(data[i], src[i]))
			return false;
	}
	return true;
}

template <typename T> const bool VectorN<T>::operator!=( const VectorN &src ) const {
	return !((*this)==src);
}


template <typename T> const bool VectorN<T>::operator==( const T d ) const {	
	for(index_type i = 0; i < _dim; ++i) {
		if(!IsAlmostSame(data[i], d))
			return false;
	}
	return true;
}

template <typename T> const bool VectorN<T>::operator!=( const T d ) const {
	return !((*this)==d);
}

//////////////////////////////////////////////////////////////////////////
// In place arithmetic
template <typename T> const VectorN<T> &VectorN<T>::operator+=(T d) { for(index_type i = 0; i < _dim; ++i) data[i] += d; return (*this); }
template <typename T> const VectorN<T> &VectorN<T>::operator-=(T d) { for(index_type i = 0; i < _dim; ++i) data[i] -= d; return (*this); }
template <typename T> const VectorN<T> &VectorN<T>::operator*=(T d) { for(index_type i = 0; i < _dim; ++i) data[i] *= d; return (*this); }
template <typename T> const VectorN<T> &VectorN<T>::operator/=(T d) { for(index_type i = 0; i < _dim; ++i) data[i] /= d; return (*this); }

//////////////////////////////////////////////////////////////////////////
// In place arithmetic; componentwise operations
template <typename T> const VectorN<T> &VectorN<T>::operator+=( const VectorN &src ) {
	__ASSERT( _dim == src.dim() );
	for(index_type i = 0; i < _dim; ++i) {
		data[i] += src[i];
	}
	return (*this);
}

template <typename T> const VectorN<T> &VectorN<T>::operator-=( const VectorN &src ) {
	__ASSERT( _dim == src.dim() );
	for(index_type i = 0; i < _dim; ++i) {
		data[i] -= src[i];
	}
	return (*this);
}

template <typename T> const VectorN<T> &VectorN<T>::operator*=( const VectorN &src ) {
	__ASSERT( _dim == src.dim() );
	for(index_type i = 0; i < _dim; ++i) {
		data[i] *= src[i];
	}
	return (*this);
}

template <typename T> const VectorN<T> &VectorN<T>::operator/=( const VectorN &src ) {
	__ASSERT( _dim == src.dim() );
	for(index_type i = 0; i < _dim; ++i) {
		data[i] /= src[i];
	}
	return (*this);
}

template <typename T> const VectorN<T> &VectorN<T>::operator+=(const Vector_V_mul_S<VectorN<T>,T> & op) {
	__ASSERT( _dim == op.V.dim() );
	for(index_type i= 0; i < _dim; ++i) {
		data[i] += (op.V[i] * op.s);
	}
	return (*this);
}	

template <typename T> const VectorN<T> &VectorN<T>::operator-=(const Vector_V_mul_S<VectorN<T>,T> & op) {
	__ASSERT( _dim == op.V.dim() );
	for(index_type i= 0; i < _dim; ++i) {
		data[i] -= (op.V[i] * op.s);
	}
	return (*this);
}

template <typename T> const VectorN<T> &VectorN<T>::operator*=(const Vector_V_mul_S<VectorN<T>,T> & op) {
	__ASSERT( _dim == op.V.dim() );
	for(index_type i= 0; i < _dim; ++i) {
		data[i] *= (op.V[i] * op.s);
	}
	return (*this);
}

template <typename T> const VectorN<T> &VectorN<T>::operator/=(const Vector_V_mul_S<VectorN<T>,T> & op) {
	__ASSERT( _dim == op.V.dim() );
	for(index_type i= 0; i < _dim; ++i) {
		data[i] /= (op.V[i] * op.s);
	}
	return (*this);
}

//////////////////////////////////////////////////////////////////////////
// IO
template <typename T> std::ostream &operator<<(std::ostream &os, const VectorN<T> &v) {
	os << "VectorN Dimension : " << v.dim() << std::endl;
	os << "[";
	for(VectorN<T>::index_type i = 0; i < v.dim() -1 ; ++i) {
		os << std::setw(8) << v[i] << ", "; 
	}
	os << std::setw(8) << v[v.dim()-1] << "]";
	return os;	
}

//////////////////////////////////////////////////////////////////////////
// Friends Operators
template <typename T> void neg( VectorN<T> &ret, const VectorN<T> &a )						{ VectorOperation<VectorN<T>>().neg(ret, a);		}
template <typename T> void add( VectorN<T> &ret, const VectorN<T> &a, const VectorN<T> &b ) { VectorOperation<VectorN<T>>().add(ret, a, b);		}
template <typename T> void dif( VectorN<T> &ret, const VectorN<T> &a, const VectorN<T> &b ) { VectorOperation<VectorN<T>>().dif(ret, a, b);		}
template <typename T> void mul( VectorN<T> &ret, const VectorN<T> &v, const T c )			{ VectorOperation<VectorN<T>>().mul(ret, v, c);		}
template <typename T> void mul( VectorN<T> &ret, const T c, const VectorN<T> &v )			{ mul(ret, v, c);								}
template <typename T> void div( VectorN<T> &ret, const VectorN<T> &v, const T c )			{ __ASSERT(c != 0); mul(ret, v, (T)(1.0f/c));	}

template <typename T> void unitVector(VectorN<T> &ret, const VectorN<T> &a) {
	T _mag = mag(a);
	if(!IsAlmostZero(_mag)) ret = a / _mag;
	else					ret.zero();
}

template <typename T> const T normalize(VectorN<T> &a) {
	T _mag = mag(a);
	if(!IsAlmostZero(_mag)) a /= _mag;
	else					a.zero();

	return _mag;
}

template <typename T> const T l1Norm(const VectorN<T> &a)	{ return VectorNorm_L1<VectorN<T>>()(a);	}
template <typename T> const T l2Norm(const VectorN<T> &a)	{ return VectorNorm_L2<VectorN<T>>()(a);	}
template <typename T> const T lInfNorm(const VectorN<T> &a) { return VectorNorm_Inf<VectorN<T>>()(a);	}

template <typename T> const T mag( const VectorN<T> &a)		{ return l2Norm(a);							}
template <typename T> const T sqrmag( const VectorN<T> &a)	{ return SQ(mag(a));						}	
template <typename T> const T sum( const VectorN<T> &a)		{ return VectorOperation::sum(a);	}

template <typename T> const T dot( const VectorN<T> &a, const VectorN<T> &b ) { __ASSERT( a.dim() == b.dim() );	return VectorOperation<VectorN<T>>().dot(a,b); }

template <typename T> const T length( const VectorN<T> &a, const VectorN<T> &b )	{ return VectorOperation<VectorN<T>>::length(a,b);		}
template <typename T> const T sqrlength( const VectorN<T> &a, const VectorN<T> &b ) { return VectorOperation<VectorN<T>>::sqrlength(a,b);	}

template <typename T> void maxVec(VectorN<T> &ret,const VectorN<T> &a,const VectorN<T> &b) { return VectorOperation<VectorN<T>>::maxVec(ret,a,b); }
template <typename T> void minVec(VectorN<T> &ret,const VectorN<T> &a,const VectorN<T> &b) { return VectorOperation<VectorN<T>>::minVec(ret,a,b); }

template <typename T> const Vector_neg_V<VectorN<T>> operator-(const VectorN<T> &a)						   { return Vector_neg_V<VectorN<T>>(a);			}
template <typename T> const Vector_V_add_V<VectorN<T>> operator+(const VectorN<T> &a, const VectorN<T> &b) { return Vector_V_add_V<VectorN<T>>(a,b);		}
template <typename T> const Vector_V_dif_V<VectorN<T>> operator-(const VectorN<T> &a, const VectorN<T> &b) { return Vector_V_dif_V<VectorN<T>>(a,b);		}
template <typename T> const Vector_V_mul_S<VectorN<T>,T> operator*(const VectorN<T> &a, const T &c)		   { return Vector_V_mul_S<VectorN<T>,T>(a,c);		}
template <typename T> const Vector_V_mul_S<VectorN<T>,T> operator*(const T &c, const VectorN<T> &a)		   { return Vector_V_mul_S<VectorN<T>,T>(a,c);		}
template <typename T> const Vector_V_mul_S<VectorN<T>,T> operator/(const VectorN<T> &a, const T &c)		   { return Vector_V_mul_S<VectorN<T>,T>(a,T(1)/c); }

typedef VectorN<char>			VectorNc;
typedef VectorN<int>			VectorNi;
typedef VectorN<float>			VectorNf;
typedef VectorN<double>			VectorNd;

#undef __ASSERT

}; // namespace VectorN

#endif
