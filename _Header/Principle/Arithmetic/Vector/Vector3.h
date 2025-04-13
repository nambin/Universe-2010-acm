#ifndef _ARITHMETIC_VECTOR3_H_
#define _ARITHMETIC_VECTOR3_H_

#include "../../../Basis/Basis/BasicLibrary.h"
#include "../../../Basis/Random/RandomGen.h"

namespace arithmetic {

#ifdef _DEBUG
#define __ASSERT(c) assert(c)
#else
#define __ASSERT(c) {}
#endif

//////////////////////////////////////////////////////////////////////////
//	Vector3 class definition
//////////////////////////////////////////////////////////////////////////
template <typename T> class Matrix3;

template <typename T> class Vector3 {
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

	static const unsigned int getDimension() { return 3; }

	//////////////////////////////////////////////////////////////////////////
	// Constructors	
	inline explicit Vector3(T d = 0);
	inline Vector3(T d0,T d1,T d2);
	inline Vector3(const Vector3&);

	template<typename S> inline explicit Vector3(const Vector3<S>&);
	template<typename S> inline explicit Vector3(const S*);  

	explicit Vector3(const basis::RandomGen<T>&);

	//////////////////////////////////////////////////////////////////////////
	// Index operators
	inline T&		operator[](index_type i)      ;		
	inline const T  operator[](index_type i) const;		

	inline T&		operator()(index_type i)      ;		
	inline const T  operator()(index_type i) const;		
	
	inline T&		x()			{ return data[0]; }
	inline T&		y()			{ return data[1]; }
	inline T&		z()			{ return data[2]; }
	inline const T	x() const	{ return data[0]; }
	inline const T	y() const	{ return data[1]; }
	inline const T	z() const	{ return data[2]; }	

	//////////////////////////////////////////////////////////////////////////
	// Assignment and set
	inline void zero();

	inline void set(T d);
	inline void set(T d0, T d1, T d2);

	template<typename S> inline void set(const Vector3<S>	&da);
	template<typename S> inline void set(const S			*da);  	

	template<typename S> inline Vector3 &operator=(const Vector3<S> &da);
	template<typename S> inline Vector3 &operator=(const S			*da);

	//////////////////////////////////////////////////////////////////////////
	// Manipulation
	inline Vector3&	normalize();
	inline const T	mag() const;
	inline const T	sqrmag() const;

	inline const bool isUnitVector() const;

	inline Matrix3<T> getDualMatrix() const;				// Implemented in Matrix3.h
	inline Matrix3<T> getTensorProduct() const;				// Implemented in Matrix3.h
	inline void getTensorProduct(Matrix3<T>&) const;		// Implemented in Matrix3.h

	//////////////////////////////////////////////////////////////////////////
	// Get Pointer
	inline const T* const getptr() const { return data; }
	   
	//////////////////////////////////////////////////////////////////////////
	// Comparison operators
	const bool operator==(const Vector3 &da) const;
	const bool operator!=(const Vector3 &da) const;

	const bool operator==(T d) const;
	const bool operator!=(T d) const;

	//////////////////////////////////////////////////////////////////////////
	// In place arithmetic
	inline Vector3 &operator+=(T d);
	inline Vector3 &operator-=(T d);
	inline Vector3 &operator*=(T d);
	inline Vector3 &operator/=(T d);

	inline Vector3 &operator+=(const Vector3 &da);
	inline Vector3 &operator-=(const Vector3 &da);
	inline Vector3 &operator*=(const Vector3 &da);
	inline Vector3 &operator/=(const Vector3 &da);

	//////////////////////////////////////////////////////////////////////////
	// Setter as maximum and minimum
	Vector3 &maxSet(const Vector3 &da);
	Vector3 &minSet(const Vector3 &da);	

	//////////////////////////////////////////////////////////////////////////
	// Component Operators
	const T minComponent() const;
	const T maxComponent() const;
	const T minAbsComponent() const;
	const T maxAbsComponent() const;
	const index_type indexOfMinComponent() const;
	const index_type indexOfMaxComponent() const;
	const index_type indexOfMinAbsComponent() const;
	const index_type indexOfMaxAbsComponent() const;

	//////////////////////////////////////////////////////////////////////////
	// IO
	template <typename T> friend std::istream& operator>>(std::istream &is,Vector3<T> &v);
	template <typename T> friend std::ostream& operator<<(std::ostream &os,const Vector3<T> &v);

	//////////////////////////////////////////////////////////////////////////
	// Friends Operators
	template <typename T> friend const Vector3<T> operator+(const Vector3<T> &a);
	template <typename T> friend const Vector3<T> operator-(const Vector3<T> &a);
	template <typename T> friend const Vector3<T> operator+(const Vector3<T> &a,const Vector3<T> &b);
	template <typename T> friend const Vector3<T> operator-(const Vector3<T> &a,const Vector3<T> &b);
	template <typename T> friend const Vector3<T> operator*(const Vector3<T> &a,const Vector3<T> &b);
	template <typename T> friend const Vector3<T> operator/(const Vector3<T> &a,const Vector3<T> &b);

	template <typename T> friend const Vector3<T> operator+(const Vector3<T> &a,T b);
	template <typename T> friend const Vector3<T> operator-(const Vector3<T> &a,T b);
	template <typename T> friend const Vector3<T> operator*(const Vector3<T> &a,T b);
	template <typename T> friend const Vector3<T> operator/(const Vector3<T> &a,T b);

	template <typename T> friend const Vector3<T> operator+(T a,const Vector3<T> &b);
	template <typename T> friend const Vector3<T> operator-(T a,const Vector3<T> &b);
	template <typename T> friend const Vector3<T> operator*(T a,const Vector3<T> &b);
	template <typename T> friend const Vector3<T> operator/(T a,const Vector3<T> &b);

	template <typename T> friend const Vector3<T> abs(const Vector3<T> &a);
	template <typename T> friend const Vector3<T> unitVector(const Vector3<T> &a);
	template <typename T> friend const T normalize(Vector3<T> &a);

	template <typename T> friend const T l1Norm(const Vector3<T> &a);
	template <typename T> friend const T l2Norm(const Vector3<T> &a);
	template <typename T> friend const T lInfNorm(const Vector3<T> &a);
	template <typename T> friend const T mag(const Vector3<T> &a);
	template <typename T> friend const T sqrmag(const Vector3<T> &a);
	template <typename T> friend const T sum(const Vector3<T> &a);	

	template <typename T> friend const index_type dominantAxis(const Vector3<T> &v);
	template <typename T> friend const index_type subinantAxis(const Vector3<T> &v);

	template <typename T> friend const T angle(const Vector3<T> &a,const Vector3<T> &b);
	template <typename T> friend const T dot(const Vector3<T> &a,const Vector3<T> &b);
	template <typename T> friend const Vector3<T> cross(const Vector3<T> &a,const Vector3<T> &b);
	template <typename T> friend const T tripleProduct(const Vector3<T> &a,const Vector3<T> &b,const Vector3<T> &c);
	template <typename T> friend const T length(const Vector3<T> &a, const Vector3<T> &b);
	template <typename T> friend const T sqrlength(const Vector3<T> &a, const Vector3<T> &b);

	template <typename T> friend const Vector3<T> reflect(const Vector3<T> &a, const Vector3<T> &b);

	template <typename T> friend const Vector3<T> maxVec(const Vector3<T> &a,const Vector3<T> &b);
	template <typename T> friend const Vector3<T> minVec(const Vector3<T> &a,const Vector3<T> &b);

public : 	
	T data[3];
};

//////////////////////////////////////////////////////////////////////////
// Constructors
template <typename T> inline Vector3<T>::Vector3(T d)					{ data[0] = data[1] = data[2] = d;					}
template <typename T> inline Vector3<T>::Vector3(T d0,T d1,T d2)		{ data[0] = d0;   data[1] = d1;   data[2] = d2;		}
template <typename T> inline Vector3<T>::Vector3(const Vector3& v)		{ data[0] = v[0]; data[1] = v[1]; data[2] = v[2];	}

template <typename T> template<typename S> inline Vector3<T>::Vector3( const Vector3<S> &da ) {
	data[0] = static_cast<T>(da[0]);
	data[1] = static_cast<T>(da[1]);
	data[2] = static_cast<T>(da[2]);
}

template <typename T> template<typename S> inline Vector3<T>::Vector3(const S *da) {
	data[0] = static_cast<T>(da[0]);
	data[1] = static_cast<T>(da[1]);
	data[2] = static_cast<T>(da[2]);
}

template <typename T> Vector3<T>::Vector3(const basis::RandomGen<T>& random) {
	data[0] = random(); 
	data[1] = random(); 
	data[2] = random();
}

//////////////////////////////////////////////////////////////////////////
// Index operators
template <typename T> inline T &Vector3<T>::operator[](index_type i)			{ __ASSERT(i<3); return data[i]; }
template <typename T> inline T &Vector3<T>::operator()(index_type i)			{ __ASSERT(i<3); return data[i]; }
template <typename T> inline const T Vector3<T>::operator[](index_type i) const { __ASSERT(i<3); return data[i]; }
template <typename T> inline const T Vector3<T>::operator()(index_type i) const { __ASSERT(i<3); return data[i]; }
  
//////////////////////////////////////////////////////////////////////////
// Assignment and set
template <typename T> inline void Vector3<T>::zero()				{ data[0] = data[1] = data[2] = 0;			}
template <typename T> inline void Vector3<T>::set(T d)				{ data[0] = data[1] = data[2] = d;			}
template <typename T> inline void Vector3<T>::set(T d0, T d1, T d2) { data[0] = d0; data[1] = d1; data[2] = d2;	}

template <typename T> template<typename S> inline void Vector3<T>::set(const Vector3<S> &da) {
	data[0] = static_cast<T>(da[0]);
	data[1] = static_cast<T>(da[1]);
	data[2] = static_cast<T>(da[2]);
}
  
template <typename T> template<typename S> inline void Vector3<T>::set(const S *da) {
	data[0] = static_cast<T>(da[0]);
	data[1] = static_cast<T>(da[1]);
	data[2] = static_cast<T>(da[2]);
}  

template <typename T> template<typename S> inline Vector3<T> &Vector3<T>::operator=(const Vector3<S> &da) { 
	data[0] = static_cast<T>(da[0]);
	data[1] = static_cast<T>(da[1]);
	data[2] = static_cast<T>(da[2]);
	return (*this); 
}

template <typename T> template<typename S> inline Vector3<T> &Vector3<T>::operator=(const S *da) { 
	data[0] = static_cast<T>(da[0]);
	data[1] = static_cast<T>(da[1]);
	data[2] = static_cast<T>(da[2]);
	return (*this); 
}

//////////////////////////////////////////////////////////////////////////
// Manipulation
template <typename T> inline Vector3<T> &Vector3<T>::normalize()	{ arithmetic::normalize(*this); return (*this); }
template <typename T> inline const T Vector3<T>::mag() const		{ return arithmetic::mag(*this);				}
template <typename T> inline const T Vector3<T>::sqrmag() const		{ return arithmetic::sqrmag(*this);				}

template <typename T> inline const bool Vector3<T>::isUnitVector() const {
	return IsAlmostSame(sqrmag(), T(1));
}
	 
//////////////////////////////////////////////////////////////////////////
// Comparison operators
template <typename T> const bool Vector3<T>::operator==(const Vector3 &da) const {
	return (basis::IsAlmostSame(data[0], da[0]) && basis::IsAlmostSame(data[1], da[1]) && basis::IsAlmostSame(data[2], da[2]));
}

template <typename T> const bool Vector3<T>::operator!=(const Vector3 &da) const {
	return !(*this == da);	
}
  
template <typename T> const bool Vector3<T>::operator==(T d) const {
	return (basis::IsAlmostSame(data[0], d) && basis::IsAlmostSame(data[1], d) && basis::IsAlmostSame(data[2], d));
}

template <typename T> const bool Vector3<T>::operator!=(T d) const {
	return !(*this == d);
}

//////////////////////////////////////////////////////////////////////////
// In place arithmetic
template <typename T> inline Vector3<T> &Vector3<T>::operator+=(T d) { data[0] += d; data[1] += d; data[2] += d; return (*this); }
template <typename T> inline Vector3<T> &Vector3<T>::operator-=(T d) { data[0] -= d; data[1] -= d; data[2] -= d; return (*this); }
template <typename T> inline Vector3<T> &Vector3<T>::operator*=(T d) { data[0] *= d; data[1] *= d; data[2] *= d; return (*this); }
template <typename T> inline Vector3<T> &Vector3<T>::operator/=(T d) { data[0] /= d; data[1] /= d; data[2] /= d; return (*this); }

//////////////////////////////////////////////////////////////////////////
// In place arithmetic; componentwise operations
template <typename T> inline Vector3<T> &Vector3<T>::operator+=(const Vector3 &da) { data[0] += da[0]; data[1] += da[1]; data[2] += da[2]; return (*this); }
template <typename T> inline Vector3<T> &Vector3<T>::operator-=(const Vector3 &da) { data[0] -= da[0]; data[1] -= da[1]; data[2] -= da[2]; return (*this); }
template <typename T> inline Vector3<T> &Vector3<T>::operator*=(const Vector3 &da) { data[0] *= da[0]; data[1] *= da[1]; data[2] *= da[2]; return (*this); }
template <typename T> inline Vector3<T> &Vector3<T>::operator/=(const Vector3 &da) { data[0] /= da[0]; data[1] /= da[1]; data[2] /= da[2]; return (*this); }

//////////////////////////////////////////////////////////////////////////
// Setter as maximum and minimum
template <typename T> Vector3<T> &Vector3<T>::maxSet(const Vector3 &da) { return (*this) = maxVec(*this, da); }
template <typename T> Vector3<T> &Vector3<T>::minSet(const Vector3 &da) { return (*this) = minVec(*this, da); }

//////////////////////////////////////////////////////////////////////////
// Component Operator
template <typename T> const T Vector3<T>::minComponent() const		{ return data[indexOfMinComponent()];	 }
template <typename T> const T Vector3<T>::maxComponent() const		{ return data[indexOfMaxComponent()];	 }
template <typename T> const T Vector3<T>::minAbsComponent() const	{ return data[indexOfMinAbsComponent()]; }
template <typename T> const T Vector3<T>::maxAbsComponent() const	{ return data[indexOfMaxAbsComponent()]; }

template <typename T> const typename Vector3<T>::index_type Vector3<T>::indexOfMinComponent() const	{
	index_type index = 0;
	T tmp = data[0];

	if(data[1] < tmp)	{	tmp = data[1];	index = 1;	}
	if(data[2] < tmp)	{	index = 2;					}

	return index;
}

template <typename T> const typename Vector3<T>::index_type Vector3<T>::indexOfMaxComponent() const		{ return (-(*this)).indexOfMinComponent();		}
template <typename T> const typename Vector3<T>::index_type Vector3<T>::indexOfMinAbsComponent() const	{ return (abs(*this)).indexOfMinComponent();	}
template <typename T> const typename Vector3<T>::index_type Vector3<T>::indexOfMaxAbsComponent() const	{ return (-abs(*this)).indexOfMinComponent();	}

//////////////////////////////////////////////////////////////////////////
// Stream
template <typename T> std::istream &operator>>(std::istream &is,Vector3<T> &v)	{
	T tmp;
	is >> tmp;	v.data[0] = tmp;
	is >> tmp;	v.data[1] = tmp;
	is >> tmp;	v.data[2] = tmp;
	return is;
}

template <typename T> std::ostream &operator<<(std::ostream &os,const Vector3<T> &v) {
	os << "[" << v[0] << ", " << v[1] << ", " << v[2] << "]";
	return os;
}

//////////////////////////////////////////////////////////////////////////
// Friends Operators
//////////////////////////////////////////////////////////////////////////

template <typename T> inline const Vector3<T> operator+(const Vector3<T> &a) { return a;							 }
template <typename T> inline const Vector3<T> operator-(const Vector3<T> &a) { return Vector3<T>(-a[0],-a[1],-a[2]); }

//////////////////////////////////////////////////////////////////////////

template <typename T> inline const Vector3<T> operator+(const Vector3<T> &a,const Vector3<T> &b) { return Vector3<T>(a[0] + b[0], a[1] + b[1], a[2] + b[2]); }
template <typename T> inline const Vector3<T> operator-(const Vector3<T> &a,const Vector3<T> &b) { return Vector3<T>(a[0] - b[0], a[1] - b[1], a[2] - b[2]); }
template <typename T> inline const Vector3<T> operator*(const Vector3<T> &a,const Vector3<T> &b) { return Vector3<T>(a[0] * b[0], a[1] * b[1], a[2] * b[2]); }
template <typename T> inline const Vector3<T> operator/(const Vector3<T> &a,const Vector3<T> &b) { return Vector3<T>(a[0] / b[0], a[1] / b[1], a[2] / b[2]); }

template <typename T> inline const Vector3<T> operator+(const Vector3<T> &a,T b) { return Vector3<T>(a[0] + b, a[1] + b, a[2] + b); }
template <typename T> inline const Vector3<T> operator-(const Vector3<T> &a,T b) { return Vector3<T>(a[0] - b, a[1] - b, a[2] - b); }
template <typename T> inline const Vector3<T> operator*(const Vector3<T> &a,T b) { return Vector3<T>(a[0] * b, a[1] * b, a[2] * b); }
template <typename T> inline const Vector3<T> operator/(const Vector3<T> &a,T b) { return Vector3<T>(a[0] / b, a[1] / b, a[2] / b); }

template <typename T> inline const Vector3<T> operator+(T a,const Vector3<T> &b) { return Vector3<T>(a + b[0], a + b[1], a + b[2]); }
template <typename T> inline const Vector3<T> operator-(T a,const Vector3<T> &b) { return Vector3<T>(a - b[0], a - b[1], a - b[2]); }
template <typename T> inline const Vector3<T> operator*(T a,const Vector3<T> &b) { return Vector3<T>(a * b[0], a * b[1], a * b[2]); }
template <typename T> inline const Vector3<T> operator/(T a,const Vector3<T> &b) { return Vector3<T>(a / b[0], a / b[1], a / b[2]); }

//////////////////////////////////////////////////////////////////////////

template <typename T> const Vector3<T> abs(const Vector3<T> &a) {
	return  Vector3<T>(((a[0]>0)?a[0]:-a[0]), ((a[1]>0)?a[1]:-a[1]),((a[2]>0)?a[2]:-a[2]));
} 

template <typename T> const Vector3<T> unitVector(const Vector3<T> &a)	{
	T m = mag(a);
	if(!basis::IsAlmostZero(m))	return a / m;
	else						return Vector3<T>(0);
}

template <typename T> const T normalize(Vector3<T> &a) {
	T m = mag(a);
	if(!basis::IsAlmostZero(m))	a /= m;
	else						a = Vector3<T>(0);

	return m;
}

//////////////////////////////////////////////////////////////////////////
// Norm
template <typename T> const T l1Norm(const Vector3<T> &a)	{ return sum(abs(a));							}
template <typename T> const T l2Norm(const Vector3<T> &a)	{ return a.mag();								}
template <typename T> const T lInfNorm(const Vector3<T> &a) { return MAX3(abs(a[0]), abs(a[1]), abs(a[2]));	}

template <typename T> const T mag(const Vector3<T> &a)				{ return static_cast<T>(sqrt(sqrmag(a)));			}
template <typename T> inline const T sqrmag(const Vector3<T> &a)	{ return (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);	}
template <typename T> inline const T sum(const Vector3<T> &a)		{ return a[0]+a[1]+a[2];							}

//////////////////////////////////////////////////////////////////////////

template <typename T> const typename Vector3<T>::index_type dominantAxis(const Vector3<T> &v) { return v.indexOfMaxAbsComponent(); }
template <typename T> const typename Vector3<T>::index_type subinantAxis(const Vector3<T> &v) { return v.indexOfMinAbsComponent(); }

//////////////////////////////////////////////////////////////////////////

template <typename T> const T angle(const Vector3<T> &a, const Vector3<T> &b)	{
	const T dotProduct = dot(a,b);
	const T AB = static_cast<T>(sqrt(dot(a,a)*dot(b,b)));
	const T angle = static_cast<T>(acos(dotProduct/AB));

	if(IsAlmostSame(dotProduct/AB, T(1)))
		return 0;

	if(_isnan(angle)) {
		std::cout << "In angle(Vector3&, Vector3&) : " << angle << std::endl;
		std::cout << "Cos : " << dotProduct/AB << std::endl;
		std::cout << "Vector : " << a << std::endl;
		std::cout << "Vector : " << b << std::endl;
		return 0;
	}

	return angle;
}

template <typename T> inline const T dot(const Vector3<T> &a,const Vector3<T> &b) {
  return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

// Cross-product of two vectors.
template <typename T> inline const Vector3<T> cross(const Vector3<T> &a,const Vector3<T> &b) {
  return Vector3<T>(a[1] * b[2] - b[1] * a[2],
		            a[2] * b[0] - b[2] * a[0],
		            a[0] * b[1] - b[0] * a[1]);
}

template <typename T> inline const T tripleProduct(const Vector3<T> &a,const Vector3<T> &b,const Vector3<T> &c) {
	return dot(cross(a,b),c);
}

template <typename T> const T length(const Vector3<T> &a, const Vector3<T> &b)	{
	return static_cast<T>(sqrt(sqrlength(a,b)));
}

template <typename T> inline const T sqrlength(const Vector3<T> &a, const Vector3<T> &b)	{
	return ((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]));
}

// b should be unit Vector
template <typename T> inline const Vector3<T> reflect(const Vector3<T> &a, const Vector3<T> &b)	{
	return (2 * dot(a,b)) * b - a;
}

//////////////////////////////////////////////////////////////////////////

template <typename T> const Vector3<T> maxVec(const Vector3<T> &a,const Vector3<T> &b) {
  return Vector3<T>((a[0]>b[0])?a[0]:b[0],
		            (a[1]>b[1])?a[1]:b[1],
		            (a[2]>b[2])?a[2]:b[2]);
}
 
template <typename T> const Vector3<T> minVec(const Vector3<T> &a,const Vector3<T> &b) {
  return Vector3<T>((a[0]<b[0])?a[0]:b[0],
		            (a[1]<b[1])?a[1]:b[1],
		            (a[2]<b[2])?a[2]:b[2]);
}

typedef Vector3<char>			Vector3c;
typedef Vector3<int>			Vector3i;
typedef Vector3<float>			Vector3f;
typedef Vector3<double>			Vector3d;

#undef __ASSERT

}; // namespace arithmetic

#endif

