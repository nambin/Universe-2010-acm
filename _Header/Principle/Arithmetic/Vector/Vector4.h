#ifndef _ARITHMETIC_VECTOR4_H_
#define _ARITHMETIC_VECTOR4_H_

#include "../../../Basis/Basis/BasicLibrary.h"
#include "Vector3.h"

namespace arithmetic {

#ifdef _DEBUG
#define __ASSERT(c) assert(c)
#else
#define __ASSERT(c) {}
#endif

//////////////////////////////////////////////////////////////////////////
//	Vector4 class definition
//////////////////////////////////////////////////////////////////////////
template <typename T> class Vector4 {
public :
	typedef unsigned int		index_type;
	typedef index_type*			index_pointer ;
	typedef const index_type*	index_const_pointer ;
	typedef index_type&			index_reference ;
	typedef const index_type&	index_const_reference ;

	typedef T					value_type;
	typedef value_type*			value_pointer ;       
	typedef const value_type*	value_const_pointer ; 
	typedef value_type&			value_reference ;     
	typedef const value_type&	value_const_reference;

	static const unsigned int getDimension() { return 4; }

	//////////////////////////////////////////////////////////////////////////
	// Constructors
	inline explicit Vector4(T d = 0);
	inline Vector4(T d0,T d1,T d2,T d3);

	inline Vector4(const Vector4 &da);
	inline Vector4(const T       *da);  
	inline Vector4(const Vector3<T> &v, T w = 1);

	//////////////////////////////////////////////////////////////////////////
	// Index operators
	inline T&		operator[](index_type i)      ;
	inline const T	operator[](index_type i) const;

	inline T&		operator()(index_type i)      ;
	inline const T  operator()(index_type i) const;

	inline const T x() const { return data[0]; }
	inline const T y() const { return data[1]; }
	inline const T z() const { return data[2]; }
	inline const T w() const { return data[3]; }
	inline const Vector3<T> xyz() const { return Vector3<T>(data[0], data[1], data[2]);	}

	//////////////////////////////////////////////////////////////////////////
	// Assignment and set
	inline void set(T d);
	inline void set(T d0, T d1, T d2, T d3);

	inline void set(const Vector4 &da);
	inline void set(const T       *da);  
		
	inline Vector4 &operator=(const Vector4 &da);
	inline Vector4 &operator=(const T       *da);  

	//////////////////////////////////////////////////////////////////////////
	// Get Pointer
	inline const T* const getptr() const { return data; }

	//////////////////////////////////////////////////////////////////////////
	// Comparison operators
	inline const bool operator==(const Vector4 &da) const;
	inline const bool operator!=(const Vector4 &da) const;

	inline const bool operator==(T d) const;
	inline const bool operator!=(T d) const;

	//////////////////////////////////////////////////////////////////////////
	// In place arithmetic
	inline Vector4 &operator+=(T d);
	inline Vector4 &operator-=(T d);
	inline Vector4 &operator*=(T d);
	inline Vector4 &operator/=(T d);

	inline Vector4 &operator+=(const Vector4 &da);
	inline Vector4 &operator-=(const Vector4 &da);
	inline Vector4 &operator*=(const Vector4 &da);
	inline Vector4 &operator/=(const Vector4 &da);

	inline const Vector3<T> getVec3() const;

	//////////////////////////////////////////////////////////////////////////
	// Friends Operators
	template <typename T> friend const Vector4<T> operator+(const Vector4<T> &a);
	template <typename T> friend const Vector4<T> operator-(const Vector4<T> &a);
	template <typename T> friend const Vector4<T> operator+(const Vector4<T> &a,const Vector4<T> &b);
	template <typename T> friend const Vector4<T> operator-(const Vector4<T> &a,const Vector4<T> &b);
	template <typename T> friend const Vector4<T> operator*(const Vector4<T> &a,const Vector4<T> &b);
	template <typename T> friend const Vector4<T> operator/(const Vector4<T> &a,const Vector4<T> &b);

	template <typename T> friend const Vector4<T> operator+(const Vector4<T> &a,T b);
	template <typename T> friend const Vector4<T> operator-(const Vector4<T> &a,T b);
	template <typename T> friend const Vector4<T> operator*(const Vector4<T> &a,T b);
	template <typename T> friend const Vector4<T> operator/(const Vector4<T> &a,T b);

	template <typename T> friend const Vector4<T> operator+(T a,const Vector4<T> &b);
	template <typename T> friend const Vector4<T> operator-(T a,const Vector4<T> &b);
	template <typename T> friend const Vector4<T> operator*(T a,const Vector4<T> &b);
	template <typename T> friend const Vector4<T> operator/(T a,const Vector4<T> &b);

	template <typename T> friend std::istream& operator>>(std::istream &is,Vector4<T> &v);
	template <typename T> friend std::ostream& operator<<(std::ostream &os,const Vector4<T> &v);

public : 	
	T data[4];
};

// Constructors
template <typename T> inline Vector4<T>::Vector4(T d)				  { data[0] = data[1] = data[2] = data[3] = d;							}
template <typename T> inline Vector4<T>::Vector4(T d0,T d1,T d2,T d3) {	data[0] = d0;	 data[1] = d1;	  data[2] = d2;	   data[3] = d3;	}
template <typename T> inline Vector4<T>::Vector4( const Vector4 &da ) {	data[0] = da[0]; data[1] = da[1]; data[2] = da[2]; data[3] = da[3]; }
template <typename T> inline Vector4<T>::Vector4(const T *da)		  {	data[0] = da[0]; data[1] = da[1]; data[2] = da[2]; data[3] = da[3]; }

template <typename T> inline Vector4<T>::Vector4( const Vector3<T> &v, T w ) { data[0] = v[0]; data[1] = v[1]; data[2] = v[2]; data[3] = w; }

// Index operators
template <typename T> inline T &Vector4<T>::operator[](index_type i)			{ __ASSERT(i<4); return data[i]; }
template <typename T> inline T &Vector4<T>::operator()(index_type i)			{ __ASSERT(i<4); return data[i]; }
template <typename T> inline const T Vector4<T>::operator[](index_type i) const { __ASSERT(i<4); return data[i]; }
template <typename T> inline const T Vector4<T>::operator()(index_type i) const { __ASSERT(i<4); return data[i]; }

// Assignment and set
template <typename T> inline void Vector4<T>::set(T d)							  { data[0] = d; data[1] = d; data[2] = d; data[3] = d;									}
template <typename T> inline void Vector4<T>::set(T d0, T d1, T d2,T d3)		  { data[0] = d0;data[1] = d1;data[2] = d2;data[3] = d3;								}
template <typename T> inline void Vector4<T>::set(const Vector4 &da)			  { data[0] = da[0]; data[1] = da[1]; data[2] = da[2]; data[3] = da[3];					}
template <typename T> inline void Vector4<T>::set(const T *da)					  { data[0] = da[0]; data[1] = da[1]; data[2] = da[2]; data[3] = da[3];					}
template <typename T> inline Vector4<T> &Vector4<T>::operator=(const Vector4 &da) {	data[0] = da[0]; data[1] = da[1]; data[2] = da[2]; data[3] = da[3]; return (*this); }
template <typename T> inline Vector4<T> &Vector4<T>::operator=(const T *da)		  { data[0] = da[0]; data[1] = da[1]; data[2] = da[2]; data[3] = da[3]; return (*this); }

// Comparison operators
template <typename T> inline const bool Vector4<T>::operator==(const Vector4 &da) const {
	return (IsAlmostSame(data[0], da[0]) && IsAlmostSame(data[1], da[1]) && IsAlmostSame(data[2], da[2]) && IsAlmostSame(data[3], da[3]));
}

template <typename T> inline const bool Vector4<T>::operator!=(const Vector4 &da) const {
	return !(*this == da);	
}

template <typename T> inline const bool Vector4<T>::operator==(T d) const {
	return (IsAlmostSame(data[0], d) && IsAlmostSame(data[1], d) && IsAlmostSame(data[2], d) && IsAlmostSame(data[3], d));
}

template <typename T> inline const bool Vector4<T>::operator!=(T d) const {
	return !(*this == d);
}

// In place arithmetic
template <typename T> inline Vector4<T> &Vector4<T>::operator+=(T d) { data[0] += d; data[1] += d; data[2] += d; data[3] += d; return (*this); }
template <typename T> inline Vector4<T> &Vector4<T>::operator-=(T d) { data[0] -= d; data[1] -= d; data[2] -= d; data[3] -= d; return (*this); }
template <typename T> inline Vector4<T> &Vector4<T>::operator*=(T d) { data[0] *= d; data[1] *= d; data[2] *= d; data[3] *= d; return (*this); }
template <typename T> inline Vector4<T> &Vector4<T>::operator/=(T d) { data[0] /= d; data[1] /= d; data[2] /= d; data[3] /= d; return (*this); }

// In place arithmetic; componentwise operations
template <typename T> inline Vector4<T> &Vector4<T>::operator+=(const Vector4 &da) { data[0] += da[0]; data[1] += da[1]; data[2] += da[2]; data[3] += da[3]; return (*this); }
template <typename T> inline Vector4<T> &Vector4<T>::operator-=(const Vector4 &da) { data[0] -= da[0]; data[1] -= da[1]; data[2] -= da[2]; data[3] -= da[3]; return (*this); }
template <typename T> inline Vector4<T> &Vector4<T>::operator*=(const Vector4 &da) { data[0] *= da[0]; data[1] *= da[1]; data[2] *= da[2]; data[3] *= da[3]; return (*this); }
template <typename T> inline Vector4<T> &Vector4<T>::operator/=(const Vector4 &da) { data[0] /= da[0]; data[1] /= da[1]; data[2] /= da[2]; data[3] /= da[3]; return (*this); }

template <typename T> inline const Vector3<T> Vector4<T>::getVec3() const {
	if(data[3] == 0.0f || data[3] == 1.0f)
		return Vector3<T>(data[0], data[1], data[2]);
	else
		return Vector3<T>(data[0]/data[3], data[1]/data[3], data[2]/data[3]);
}

//////////////////////////////////////////////////////////////////////////
// Friends Operators
//////////////////////////////////////////////////////////////////////////

template <typename T> inline const Vector4<T> operator+(const Vector4<T> &a) { return a;									}
template <typename T> inline const Vector4<T> operator-(const Vector4<T> &a) { return Vector4<T>(-a[0],-a[1],-a[2],-a[3]);  }

//////////////////////////////////////////////////////////////////////////
template <typename T> inline const Vector4<T> operator+(const Vector4<T> &a,const Vector4<T> &b) { return Vector4<T>(a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]); }
template <typename T> inline const Vector4<T> operator-(const Vector4<T> &a,const Vector4<T> &b) { return Vector4<T>(a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3]); }
template <typename T> inline const Vector4<T> operator*(const Vector4<T> &a,const Vector4<T> &b) { return Vector4<T>(a[0] * b[0], a[1] * b[1], a[2] * b[2], a[3] * b[3]); }
template <typename T> inline const Vector4<T> operator/(const Vector4<T> &a,const Vector4<T> &b) { return Vector4<T>(a[0] / b[0], a[1] / b[1], a[2] / b[2], a[3] / b[3]); }

template <typename T> inline const Vector4<T> operator+(const Vector4<T> &a,T b) { return Vector4<T>(a[0] + b, a[1] + b, a[2] + b, a[3] + b); }
template <typename T> inline const Vector4<T> operator-(const Vector4<T> &a,T b) { return Vector4<T>(a[0] - b, a[1] - b, a[2] - b, a[3] - b); }
template <typename T> inline const Vector4<T> operator*(const Vector4<T> &a,T b) { return Vector4<T>(a[0] * b, a[1] * b, a[2] * b, a[3] * b); }
template <typename T> inline const Vector4<T> operator/(const Vector4<T> &a,T b) { return Vector4<T>(a[0] / b, a[1] / b, a[2] / b, a[3] / b); }

template <typename T> inline const Vector4<T> operator+(T a,const Vector4<T> &b) { return Vector4<T>(a + b[0], a + b[1], a + b[2], a + b[3]); }
template <typename T> inline const Vector4<T> operator-(T a,const Vector4<T> &b) { return Vector4<T>(a - b[0], a - b[1], a - b[2], a - b[3]); }
template <typename T> inline const Vector4<T> operator*(T a,const Vector4<T> &b) { return Vector4<T>(a * b[0], a * b[1], a * b[2], a * b[3]); }
template <typename T> inline const Vector4<T> operator/(T a,const Vector4<T> &b) { return Vector4<T>(a / b[0], a / b[1], a / b[2], a / b[3]); }

///////////////////////////////////////////////////////////////
// Stream
template <typename T> std::istream &operator>>(std::istream &is,Vector4<T> &v)	{
	T tmp;
	is >> tmp;	v.data[0] = tmp;
	is >> tmp;	v.data[1] = tmp;
	is >> tmp;	v.data[2] = tmp;
	is >> tmp;	v.data[3] = tmp;
	return is;
}
template <typename T> std::ostream &operator<<(std::ostream &os,const Vector4<T> &v) {
	os << "[" << v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3] << "]";
	return os;
}

typedef Vector4<char>			Vector4c;
typedef Vector4<int>			Vector4i;
typedef Vector4<float>			Vector4f;
typedef Vector4<double>			Vector4d;

#undef __ASSERT

}; // namespace arithmetic

#endif
