#ifndef _ARITHMETIC_VECTOR2_H_
#define _ARITHMETIC_VECTOR2_H_

#include "../../../Basis/Basis/BasicLibrary.h"
#include "../../../Basis/Random/RandomGen.h"

namespace arithmetic {

#ifdef _DEBUG
#define __ASSERT(c) assert(c)
#else
#define __ASSERT(c) {}
#endif

//////////////////////////////////////////////////////////////////////////
//	Vector2 class definition
//////////////////////////////////////////////////////////////////////////
template <typename T> class Vector2 {
public:
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

	static const unsigned int getDimension() { return 2; }

	//////////////////////////////////////////////////////////////////////////
	// Constructors	
	inline explicit Vector2(T d = 0);	
	inline Vector2(T d0, T d1);	
	inline Vector2(const Vector2<T>&);

	template<typename S> inline explicit Vector2(const Vector2<S>&	da);
	template<typename S> inline explicit Vector2(const S*			da);  

	inline explicit Vector2(const basis::RandomGen<T>&);

	//////////////////////////////////////////////////////////////////////////
	// Index operators
	inline T&		operator[](index_type i)      ;		
	inline const T  operator[](index_type i) const;		

	inline T&		operator()(index_type i)      ;		
	inline const T  operator()(index_type i) const;		

	inline T&	 	x() 	  { return data[0]; }
	inline T&		y() 	  { return data[1]; }
	inline const T 	x() const { return data[0]; }
	inline const T 	y() const { return data[1]; }

	//////////////////////////////////////////////////////////////////////////
	// Assignment and set
	inline void zero();

	inline void set(T d);
	inline void set(T d0, T d1);				

	template<typename S> inline void set(const Vector2<S>	&da);
	template<typename S> inline void set(const S			*da);  	
					
	template<typename S> inline Vector2 &operator=(const Vector2<S> &da);
	template<typename S> inline Vector2 &operator=(const S			*da);

	//////////////////////////////////////////////////////////////////////////
	// Manipulation
	inline Vector2&	normalize();
	inline const T	mag() const;
	inline const T	sqrmag() const;

	//////////////////////////////////////////////////////////////////////////
	// Get Pointer
	inline const T* const getptr() const { return data; }

	//////////////////////////////////////////////////////////////////////////
	// Comparison operators
	const bool operator==(const Vector2 &da) const;	
	const bool operator!=(const Vector2 &da) const;	

	const bool operator==(T d) const;					
	const bool operator!=(T d) const;					

	//////////////////////////////////////////////////////////////////////////
	// In place arithmetic
	inline Vector2 &operator+=(T d);	
	inline Vector2 &operator-=(T d);	
	inline Vector2 &operator*=(T d);	
	inline Vector2 &operator/=(T d);	

	inline Vector2 &operator+=(const Vector2 &da); 
	inline Vector2 &operator-=(const Vector2 &da); 
	inline Vector2 &operator*=(const Vector2 &da); 
	inline Vector2 &operator/=(const Vector2 &da); 

	//////////////////////////////////////////////////////////////////////////
	// Setter as maximum and minimum
	Vector2 &maxSet(const Vector2 &da);	
	Vector2 &minSet(const Vector2 &da);	

	//////////////////////////////////////////////////////////////////////////
	// Static methods
	static const index_type cycleAxis(index_type axis, int direction);

	//////////////////////////////////////////////////////////////////////////
	// IO
	template <typename T> friend std::istream& operator>>(std::istream &is, Vector2<T> &v);
	template <typename T> friend std::ostream& operator<<(std::ostream &os, const Vector2<T> &v);

	//////////////////////////////////////////////////////////////////////////
	// Friends Operators
	template <typename T> friend const Vector2<T> operator+(const Vector2<T> &a);
	template <typename T> friend const Vector2<T> operator-(const Vector2<T> &a);
	template <typename T> friend const Vector2<T> operator+(const Vector2<T> &a,const Vector2<T> &b);
	template <typename T> friend const Vector2<T> operator-(const Vector2<T> &a,const Vector2<T> &b);
	template <typename T> friend const Vector2<T> operator*(const Vector2<T> &a,const Vector2<T> &b);
	template <typename T> friend const Vector2<T> operator/(const Vector2<T> &a,const Vector2<T> &b);

	template <typename T> friend const Vector2<T> operator+(const Vector2<T> &a,T b);
	template <typename T> friend const Vector2<T> operator-(const Vector2<T> &a,T b);
	template <typename T> friend const Vector2<T> operator*(const Vector2<T> &a,T b);
	template <typename T> friend const Vector2<T> operator/(const Vector2<T> &a,T b);

	template <typename T> friend const Vector2<T> operator+(T a,const Vector2<T> &b);
	template <typename T> friend const Vector2<T> operator-(T a,const Vector2<T> &b);
	template <typename T> friend const Vector2<T> operator*(T a,const Vector2<T> &b);
	template <typename T> friend const Vector2<T> operator/(T a,const Vector2<T> &b);

	template <typename T> friend const Vector2<T> abs(const Vector2<T> &a);
	template <typename T> friend const Vector2<T> unitVector(const Vector2<T> &a);
	template <typename T> friend const T normalize(const Vector2<T> &a);	

	template <typename T> friend const T l1Norm(const Vector2<T> &a);
	template <typename T> friend const T l2Norm(const Vector2<T> &a);
	template <typename T> friend const T lInfNorm(const Vector2<T> &a);
	template <typename T> friend const T mag(const Vector2<T> &a);
	template <typename T> friend const T sqrmag(const Vector2<T> &a);
	template <typename T> friend const T sum(const Vector2<T> &a);	

	template <typename T> friend const index_type dominantAxis(const Vector2<T> &v);
	template <typename T> friend const index_type subinantAxis(const Vector2<T> &v);

	template <typename T> friend const T dot(const Vector2<T> &a,const Vector2<T> &b);
	template <typename T> friend const T cross(const Vector2<T> &a,const Vector2<T> &b);
	template <typename T> friend const T length(const Vector2<T> &a, const Vector2<T> &b);
	template <typename T> friend const T sqrlength(const Vector2<T> &a, const Vector2<T> &b);

	template <typename T> friend const Vector2<T> maxVec(const Vector2<T> &a,const Vector2<T> &b);
	template <typename T> friend const Vector2<T> minVec(const Vector2<T> &a,const Vector2<T> &b);

public :
	T data[2];
};
	
//////////////////////////////////////////////////////////////////////////
// Constructors
template<typename T> inline Vector2<T>::Vector2(T d)					{ data[0] = data[1] = d;			}
template<typename T> inline Vector2<T>::Vector2(T d0,T d1)				{ data[0] = d0;	  data[1] = d1;		}
template<typename T> inline Vector2<T>::Vector2(const Vector2<T>& v)	{ data[0] = v[0]; data[1] = v[1];	}

template <typename T> template<typename S> inline Vector2<T>::Vector2( const Vector2<S> &da ) { data[0] = static_cast<T>(da[0]); data[1] = static_cast<T>(da[1]); }
template <typename T> template<typename S> inline Vector2<T>::Vector2(const S *da)			  { data[0] = static_cast<T>(da[0]); data[1] = static_cast<T>(da[1]); }

template<typename T> inline Vector2<T>::Vector2(const basis::RandomGen<T>& random) { data[0] = random(); data[1] = random(); }

//////////////////////////////////////////////////////////////////////////
// Index operators
template<typename T> inline T &Vector2<T>::operator[](index_type i)				{ __ASSERT(i<2); return data[i]; }
template<typename T> inline T &Vector2<T>::operator()(index_type i)				{ __ASSERT(i<2); return data[i]; }
template<typename T> inline const T Vector2<T>::operator[](index_type i) const	{ __ASSERT(i<2); return data[i]; }
template<typename T> inline const T Vector2<T>::operator()(index_type i) const	{ __ASSERT(i<2); return data[i]; }

//////////////////////////////////////////////////////////////////////////
// Assignment and set
template <typename T> inline void Vector2<T>::zero()			{ data[0] = data[1] = 0;		}
template <typename T> inline void Vector2<T>::set(T d)			{ data[0] = data[1] = d;		}
template <typename T> inline void Vector2<T>::set(T d0, T d1)	{ data[0] = d0;	 data[1] = d1;	}

template <typename T> template<typename S> inline void Vector2<T>::set(const Vector2<S> &da) { data[0] = static_cast<T>(da[0]);	data[1] = static_cast<T>(da[1]); }
template <typename T> template<typename S> inline void Vector2<T>::set(const S *da)			 { data[0] = static_cast<T>(da[0]);	data[1] = static_cast<T>(da[1]); }  

template <typename T> template<typename S> inline Vector2<T> &Vector2<T>::operator=(const Vector2<S> &da) { 
	data[0] = static_cast<T>(da[0]); 
	data[1] = static_cast<T>(da[1]); 
	return (*this); 
}

template <typename T> template<typename S> inline Vector2<T> &Vector2<T>::operator=(const S *da) { 
	data[0] = static_cast<T>(da[0]); 
	data[1] = static_cast<T>(da[1]); 
	return (*this); 
}

//////////////////////////////////////////////////////////////////////////
// Manipulation
template <typename T> inline Vector2<T> &Vector2<T>::normalize()	{ arithmetic::normalize(*this); return (*this); }
template <typename T> inline const T Vector2<T>::mag() const		{ return arithmetic::mag(*this);				}
template <typename T> inline const T Vector2<T>::sqrmag() const		{ return arithmetic::sqrmag(*this);				}

//////////////////////////////////////////////////////////////////////////
// Comparison operators
template <typename T> const bool Vector2<T>::operator==(const Vector2 &da) const { return (IsAlmostSame(data[0], da[0]) && IsAlmostSame(data[1],da[1]));	}
template <typename T> const bool Vector2<T>::operator!=(const Vector2 &da) const { return !(*this == da);													}
template <typename T> const bool Vector2<T>::operator==(T d) const				 { return (IsAlmostSame(data[0], d) && IsAlmostSame(data[1],d));			}
template <typename T> const bool Vector2<T>::operator!=(T d) const				 { return !(*this == d);													}

//////////////////////////////////////////////////////////////////////////
// In place arithmetic
template <typename T> inline Vector2<T> &Vector2<T>::operator+=(T d) { data[0] += d; data[1] += d; return (*this); }
template <typename T> inline Vector2<T> &Vector2<T>::operator-=(T d) { data[0] -= d; data[1] -= d; return (*this); }
template <typename T> inline Vector2<T> &Vector2<T>::operator*=(T d) { data[0] *= d; data[1] *= d; return (*this); }
template <typename T> inline Vector2<T> &Vector2<T>::operator/=(T d) { data[0] /= d; data[1] /= d; return (*this); }

//////////////////////////////////////////////////////////////////////////
// In place arithmetic; componentwise operations
template <typename T> inline Vector2<T> &Vector2<T>::operator+=(const Vector2 &da) { data[0] += da[0]; data[1] += da[1]; return (*this); }
template <typename T> inline Vector2<T> &Vector2<T>::operator-=(const Vector2 &da) { data[0] -= da[0]; data[1] -= da[1]; return (*this); }
template <typename T> inline Vector2<T> &Vector2<T>::operator*=(const Vector2 &da) { data[0] *= da[0]; data[1] *= da[1]; return (*this); }
template <typename T> inline Vector2<T> &Vector2<T>::operator/=(const Vector2 &da) { data[0] /= da[0]; data[1] /= da[1]; return (*this); }

//////////////////////////////////////////////////////////////////////////
// Setter as maximum and minimum
template <typename T> Vector2<T> &Vector2<T>::maxSet(const Vector2 &da) { return (*this) = maxVec(*this, da); }
template <typename T> Vector2<T> &Vector2<T>::minSet(const Vector2 &da) { return (*this) = minVec(*this, da); }

// Static methods
template <typename T> const typename Vector2<T>::index_type Vector2<T>::cycleAxis(index_type axis, int direction) {
	switch (axis+direction) {
		case 0: case 2: case 4: return 0;
		case 1: case 3: case 5: return 1;
		default: return (axis+direction) % 2;
	}
}

//////////////////////////////////////////////////////////////////////////
// Stream
template <typename T> std::istream &operator>>(std::istream &is, Vector2<T> &v)	{
	T tmp;
	is >> tmp;	v.data[0] = tmp;
	is >> tmp;	v.data[1] = tmp;
	return is;
}

template <typename T> std::ostream &operator<<(std::ostream &os, const Vector2<T> &v) {
	os << "[" << v[0] << ", " << v[1] << "]";
	return os;
}

//////////////////////////////////////////////////////////////////////////
// Friends Operators
//////////////////////////////////////////////////////////////////////////

template <typename T> inline const Vector2<T> operator+(const Vector2<T> &a) { return a;						}
template <typename T> inline const Vector2<T> operator-(const Vector2<T> &a) { return Vector2<T>(-a[0],-a[1]);	}

//////////////////////////////////////////////////////////////////////////

template <typename T> inline const Vector2<T> operator+(const Vector2<T> &a,const Vector2<T> &b) { return Vector2<T>(a[0] + b[0], a[1] + b[1]); }
template <typename T> inline const Vector2<T> operator-(const Vector2<T> &a,const Vector2<T> &b) { return Vector2<T>(a[0] - b[0], a[1] - b[1]); }
template <typename T> inline const Vector2<T> operator*(const Vector2<T> &a,const Vector2<T> &b) { return Vector2<T>(a[0] * b[0], a[1] * b[1]); }
template <typename T> inline const Vector2<T> operator/(const Vector2<T> &a,const Vector2<T> &b) { return Vector2<T>(a[0] / b[0], a[1] / b[1]); }

template <typename T> inline const Vector2<T> operator+(const Vector2<T> &a,T b) { return Vector2<T>(a[0] + b, a[1] + b); }
template <typename T> inline const Vector2<T> operator-(const Vector2<T> &a,T b) { return Vector2<T>(a[0] - b, a[1] - b); }
template <typename T> inline const Vector2<T> operator*(const Vector2<T> &a,T b) { return Vector2<T>(a[0] * b, a[1] * b); }
template <typename T> inline const Vector2<T> operator/(const Vector2<T> &a,T b) { return Vector2<T>(a[0] / b, a[1] / b); }

template <typename T> inline const Vector2<T> operator+(T a,const Vector2<T> &b) { return Vector2<T>(a + b[0], a + b[1]); }
template <typename T> inline const Vector2<T> operator-(T a,const Vector2<T> &b) { return Vector2<T>(a - b[0], a - b[1]); }
template <typename T> inline const Vector2<T> operator*(T a,const Vector2<T> &b) { return Vector2<T>(a * b[0], a * b[1]); }
template <typename T> inline const Vector2<T> operator/(T a,const Vector2<T> &b) { return Vector2<T>(a / b[0], a / b[1]); }

//////////////////////////////////////////////////////////////////////////

template <typename T> const Vector2<T> abs(const Vector2<T> &a) {
	return Vector2<T>(((a[0]>0)?a[0]:-a[0]), ((a[1]>0)?a[1]:-a[1]));
} 

template <typename T> const Vector2<T> unitVector(const Vector2<T> &a)	{
	T m = mag(a);
	if(!IsAlmostZero(m)) return a / m;
	else				 return Vector2<T>(0);
}

template <typename T> const T normalize(Vector2<T> &a) {
	T m = mag(a);
	if (!IsAlmostZero(m)) 
		a /= m;
	else	
		a = Vector2<T>(0);
	
	return m;
}

//////////////////////////////////////////////////////////////////////////
// Norm
template <typename T> const T l1Norm(const Vector2<T> &a)	{ return sum(abs(a));				}
template <typename T> const T l2Norm(const Vector2<T> &a)	{ return mag(a);					}
template <typename T> const T lInfNorm(const Vector2<T> &a) { return MAX(abs(a[0]), abs(a[1])); }

template <typename T> const T mag(const Vector2<T> &a)				{ return static_cast<T>(sqrt(a.sqrmag()));	}
template <typename T> inline const T sqrmag(const Vector2<T> &a)	{ return (a[0] * a[0] + a[1] * a[1]);		}
template <typename T> inline const T sum(const Vector2<T> &a)		{ return a[0]+a[1];							}

//////////////////////////////////////////////////////////////////////////

template <typename T> const typename Vector2<T>::index_type dominantAxis(const Vector2<T> &v) {
	T x,y;
	if (v[0]>0) x = v[0]; else x = -v[0];
	if (v[1]>0) y = v[1]; else y = -v[1];
	return ( x > y ) ? 0 : 1;
}

template <typename T> const typename Vector2<T>::index_type subinantAxis(const Vector2<T> &v) {
	T x,y;
	if (v[0]>0) x = v[0]; else x = -v[0];
	if (v[1]>0) y = v[1]; else y = -v[1];
	return ( x < y ) ? 0 : 1;
}

//////////////////////////////////////////////////////////////////////////
// Dot product of two vectors.
template <typename T> inline const T dot(const Vector2<T> &a,const Vector2<T> &b)	{ return (a[0] * b[0] + a[1] * b[1]);			}
template <typename T> inline const T cross(const Vector2<T> &a,const Vector2<T> &b) { return (a[0] * b[1] - b[0] * a[1]);			}

template <typename T> const T length(const Vector2<T> &a, const Vector2<T> &b)			 { return static_cast<T>(sqrt(sqrlength(a,b)));				}
template <typename T> inline const T sqrlength(const Vector2<T> &a, const Vector2<T> &b) { return (a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1]);	}

//////////////////////////////////////////////////////////////////////////
template <typename T> const Vector2<T> maxVec(const Vector2<T> &a,const Vector2<T> &b) { return Vector2<T>((a[0]>b[0])?a[0]:b[0], (a[1]>b[1])?a[1]:b[1]); }
template <typename T> const Vector2<T> minVec(const Vector2<T> &a,const Vector2<T> &b) { return Vector2<T>((a[0]<b[0])?a[0]:b[0], (a[1]<b[1])?a[1]:b[1]); }

typedef Vector2<char>			Vector2c;
typedef Vector2<int>			Vector2i;
typedef Vector2<float>			Vector2f;
typedef Vector2<double>			Vector2d;

#undef __ASSERT

};	// namespace arithmetic

#endif
