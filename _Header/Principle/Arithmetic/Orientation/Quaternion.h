#ifndef _ARITHMETIC_QUARTERNION_H_
#define _ARITHMETIC_QUARTERNION_H_

#include "../../Principle.h"
#include "../Vector/Vector3.h"
#include "../Matrix/Matrix3.h"
#include "../Matrix/Matrix4.h"

namespace arithmetic {
	
using arithmetic::Vector3;
using arithmetic::Matrix3;
using arithmetic::Matrix4;

template<typename T> class AxisAngle;
template<typename T> class EulerAngle;

//////////////////////////////////////////////////////////////////////////
// Class for Unit Quaternion (See inverse Function)
//////////////////////////////////////////////////////////////////////////
template<typename T>
class Quaternion {
public:
	//////////////////////////////////////////////////////////////////////////
	// Constructor
	Quaternion()								{ set(1.0f, Vector3<T>(0));	}
	Quaternion(T __w, T __x, T __y, T __z)		{ set(__w, __x, __y, __z);	}
	explicit Quaternion(const Vector3<T> &__v)	{ set(__v);					}	// used as Vector3
	explicit Quaternion(const Matrix3<T> &mat)	{ set(mat);					}
	explicit Quaternion(const Matrix4<T> &mat)	{ set(mat);					}
	Quaternion(const Vector3<T> &axis, T angle)	{ set(axis, angle);			}
	explicit Quaternion(const AxisAngle<T> &aa)	{ set(aa);					}
	explicit Quaternion(const EulerAngle<T> &ea){ set(ea);					}

private : 
	Quaternion(T __w, const Vector3<T> &__v)	{ set(__w, __v);		}
	void set(T __w, const Vector3<T> &__v)		{ _w = __w; _v = __v;	}

public : 
	//////////////////////////////////////////////////////////////////////////
	// Getter
	inline const T	w() const	{ return _w;	}
	inline const T	x() const	{ return _v[0];	}
	inline const T	y() const	{ return _v[1];	}
	inline const T	z() const	{ return _v[2];	}
	inline T&		w()			{ return _w;	}
	inline T&		x()			{ return _v[0];	}
	inline T&		y()			{ return _v[1];	}
	inline T&		z()			{ return _v[2];	}

	inline const Vector3<T>&	vec() const { return _v; }
	inline Vector3<T>&			vec()		{ return _v; }
	inline const Vector3<T>&	v() const	{ return _v; }
	inline Vector3<T>&			v()			{ return _v; }

	void getUnitAxisAngle(Vector3<T> &axis, T &angle);

	//////////////////////////////////////////////////////////////////////////
	// Setter
	void set(const Quaternion &qt)				{ (*this) = qt;							}
	void set(T __w, T __x, T __y, T __z)		{ set(__w, Vector3<T>(__x, __y, __z));	}
	void set(const Vector3<T> &__v)				{ _w = 0; _v = __v;						}
	void set(const Matrix3<T> &mat)				{ set(Matrix4<T>(mat));					}
	void set(const Matrix4<T> &mat);
	void set(const Vector3<T> &axis, T angle)	{ parameterize(axis, angle);			}
	inline void set(const AxisAngle<T> &aa);	// Implemented in AxisAngle.h
	inline void set(const EulerAngle<T> &ea);	// Implemented in EulerAngle.h

	inline Quaternion& setCloser(const Quaternion& qt);					

	//////////////////////////////////////////////////////////////////////////
	// In place arithmetic
	inline Quaternion& operator*=(const Quaternion &qt);
	inline Quaternion& operator/=(const Quaternion &qt);
	inline Quaternion& operator*=(T c);
	inline Quaternion& operator/=(T c);

	//////////////////////////////////////////////////////////////////////////
	// Manipulation
	Quaternion& parameterize(const Vector3<T> &axis, T theta);	

	Quaternion& normalize();
	inline Quaternion& invert();
	inline const Quaternion getInverse() const;

	const T mag() const;	
	inline const T sqrmag() const;	

	inline const bool isUnit() const;

	//////////////////////////////////////////////////////////////////////////
	// Matrix Operation
	const Matrix3<T> getMatrix3() const;
	const Matrix4<T> getMatrix4() const { return Matrix4(getMatrix3()); }											

	//////////////////////////////////////////////////////////////////////////
	// IO
	template <typename T> friend std::ostream& operator<<(std::ostream &os, const Quaternion<T> &q);

	//////////////////////////////////////////////////////////////////////////
	// Friends	
	template <typename T> friend const Quaternion<T> operator+(const Quaternion<T>& a);
	template <typename T> friend const Quaternion<T> operator-(const Quaternion<T>& a);
	template <typename T> friend const Quaternion<T> operator+(const Quaternion<T>& a, const Quaternion<T> &b);
	template <typename T> friend const Quaternion<T> operator-(const Quaternion<T>& a, const Quaternion<T> &b);
	template <typename T> friend const Quaternion<T> operator*(const Quaternion<T>& a, const Quaternion<T> &b);
	template <typename T> friend const Quaternion<T> operator/(const Quaternion<T>& a, const Quaternion<T> &b);
	template <typename T> friend const Quaternion<T> operator*(const Quaternion<T>& a, T c);
	template <typename T> friend const Quaternion<T> operator*(T c, const Quaternion<T>& a);
	template <typename T> friend const Quaternion<T> operator/(const Quaternion<T>& a, T c);

	template <typename T> friend const Quaternion<T> inverse(const Quaternion<T>& a);
	template <typename T> friend const T dot(const Quaternion<T>& a, const Quaternion<T>& b);
	template <typename T> friend const T mag(const Quaternion<T>& a);
	template <typename T> friend const T sqrmag(const Quaternion<T>& a);

	template <typename T> friend const Vector3<T> ln(const Quaternion<T> &qt);
	template <typename T> friend const Quaternion<T> exp(const Vector3<T> &v);
	template <typename T> friend const Quaternion<T> pow(const Vector3<T> &v, T c);

	template <typename T> friend const Quaternion<T> slerp(const Quaternion<T>& a, const Quaternion<T>& b, T c);
	template <typename T> friend const T distance(const Quaternion<T>& a, const Quaternion<T>& b);
	template <typename T> friend const Vector3<T> difference(const Quaternion<T>& a, const Quaternion<T>& b);

private:
	T			_w;
	Vector3<T>	_v;
};

//////////////////////////////////////////////////////////////////////////
// Getter
template<typename T>
void Quaternion<T>::getUnitAxisAngle(Vector3<T> &axis, T &angle) {
	angle = 2.0f * acos(_w);

	if(IsAlmostZero(_v.sqrmag())) axis = Vector3<T>(1,0,0);
	else						  axis = unitVector(_v);
}

//////////////////////////////////////////////////////////////////////////
// Setter
template<typename T>
void Quaternion<T>::set(const Matrix4<T> &mat) {
	Quaternion<T> &q = (*this);

	T tr, s;
	int i, j, k;
	static int next[3] = { 1, 2, 0 };

	tr = mat(0,0) + mat(1,1) + mat(2,2);
	if ( tr > 0.0f )	{
		s = sqrt( tr + 1.0f );
		q._w = ( s * 0.5f );
		s = 0.5f / s;
		q._v[0] = ( mat(1,2) - mat(2,1) ) * s;
		q._v[1] = ( mat(2,0) - mat(0,2) ) * s;
		q._v[2] = ( mat(0,1) - mat(1,0) ) * s;
	}
	else	{
		i = 0;
		if ( mat(1,1) > mat(0,0) ) i = 1;
		if ( mat(2,2) > mat(i,i) ) i = 2;

		j = next[i];
		k = next[j];

		s = sqrt( (mat(i,i)	- (mat(j,j) + mat(k,k))) + 1.0f );
		q._v[i] = s * 0.5f;
		s = 0.5f / s;
		q._w    = ( mat(j,k) - mat(k,j) ) * s;
		q._v[j] = ( mat(i,j) + mat(j,i) ) * s;
		q._v[k] = ( mat(i,k) + mat(k,i) ) * s;
	}
}

template<typename T>
inline Quaternion<T>& Quaternion<T>::setCloser(const Quaternion<T>& qt) {

	T dist_with_original = arithmetic::sqrmag((*this) - qt);		
	T dist_with_opposite = arithmetic::sqrmag((*this) + qt);

	if( dist_with_opposite < dist_with_original )	{
		_w *= -1;
		_v *= -1;
	}

	return (*this);
}

//////////////////////////////////////////////////////////////////////////
// In place arithmetic
template<typename T>
inline Quaternion<T>& Quaternion<T>::operator*=(const Quaternion<T> &qt) {
	return ((*this) = (*this) * qt);
}

template<typename T>
inline Quaternion<T>& Quaternion<T>::operator/=(const Quaternion<T> &qt) {
	return ((*this) = (*this) / qt);
}

template<typename T>
inline Quaternion<T>& Quaternion<T>::operator*=(T c) {
	return ((*this) = (*this) * c);
}

template<typename T>
inline Quaternion<T>& Quaternion<T>::operator/=(T c) {
	return ((*this) = (*this) * (1.0f/c));
}

//////////////////////////////////////////////////////////////////////////
// Manipulation
template<typename T>
Quaternion<T>& Quaternion<T>::parameterize(const Vector3<T> &axis, T theta){
	_w = cos(theta*0.5f);
	_v = unitVector(axis)*sin(theta*0.5f);
	return (*this);
}

template<typename T>
Quaternion<T>& Quaternion<T>::normalize() {
	return ((*this) /= this->mag());
}

template<typename T>
inline Quaternion<T>& Quaternion<T>::invert() {
	return ((*this) = arithmetic::inverse(*this));
}

template<typename T>
inline const Quaternion<T> Quaternion<T>::getInverse() const {
	return arithmetic::inverse(*this);
}

template<typename T>
const T Quaternion<T>::mag() const {
	return sqrt(sqrmag());
}

template<typename T>
inline const T Quaternion<T>::sqrmag() const {
	return _w*_w + _v.sqrmag();
}

template<typename T>
inline const bool Quaternion<T>::isUnit() const {
	//return IsAlmostSame((T)1.0f, sqrmag());
	return (::abs(sqrmag() - 1.0f) < 1.0e-4);
}

//////////////////////////////////////////////////////////////////////////
// Matrix Operation
template<typename T>
const Matrix3<T> Quaternion<T>::getMatrix3() const {
	Matrix3<T> m;
	T s, xs, ys, zs, wx, wy, wz, xx, xy, xz, yy, yz, zz;

	s  = 2.0f / mag();
	xs = x() * s;  ys = y() * s;  zs = z() * s;
	wx = w() * xs; wy = w() * ys; wz = w() * zs;
	xx = x() * xs; xy = x() * ys; xz = x() * zs;
	yy = y() * ys; yz = y() * zs; zz = z() * zs;

	m(0,0) = 1.0f - (yy + zz);
	m(1,0) = xy - wz;
	m(2,0) = xz + wy;
	m(0,1) = xy + wz;
	m(1,1) = 1.0f - (xx + zz);
	m(2,1) = yz - wx;
	m(0,2) = xz - wy;
	m(1,2) = yz + wx;
	m(2,2) = 1.0f - (xx + yy);

	return m;
}

//////////////////////////////////////////////////////////////////////////
// IO
template<typename T>
inline std::ostream& operator<<(std::ostream &os, const Quaternion<T> &q) {
	os << "[" << q._w << ", " << q._v << "]";
	return os;
}

//////////////////////////////////////////////////////////////////////////
// Friends Functions
template<typename T>
inline const Quaternion<T> operator+(const Quaternion<T>& a) {	
	return a;
}

template<typename T>
inline const Quaternion<T> operator-(const Quaternion<T>& a) {	
	return Quaternion<T>(-a._w, -a._v);
}

template<typename T>
inline const Quaternion<T> operator+(const Quaternion<T>& a, const Quaternion<T> &b) {	
	return Quaternion<T>(a._w+b._w, a._v+b._v);
}

template<typename T>
inline const Quaternion<T> operator-(const Quaternion<T>& a, const Quaternion<T> &b) {	
	return Quaternion<T>(a._w-b._w, a._v-b._v);
}

template<typename T>
inline const Quaternion<T> operator*(const Quaternion<T>& a, const Quaternion<T> &b) {	
	return Quaternion<T>(a._w*b._w-dot(a._v, b._v), a._w*b._v+b._w*a._v+cross(a._v,b._v));
}

template<typename T>
inline const Quaternion<T> operator/(const Quaternion<T>& a, const Quaternion<T> &b) {	
	return a * inverse(b);
}

template<typename T>
inline const Quaternion<T> operator*(const Quaternion<T>& a, T c) {	
	return Quaternion<T>(a._w*c, a._v*c);
}

template<typename T>
inline const Quaternion<T> operator*(T c, const Quaternion<T>& a) {	
	return a*c;
}

template<typename T>
inline const Quaternion<T> operator/(const Quaternion<T>& a, T c) {	
	return a * (1.0f/c);
}

template<typename T>
inline const Quaternion<T> inverse(const Quaternion<T>& qt) {		// Assume that Quaternion<T> is Unit Quaternion
	return Quaternion<T>(qt._w, -qt._v);
}

template<typename T>
inline const T dot(const Quaternion<T>& a, const Quaternion<T>& b) {	
	return a._w * b._w + dot(a._v, b._v);
}

template<typename T>
inline const T mag(const Quaternion<T>& a) {
	return a.mag();
}

template<typename T>
inline const T sqrmag(const Quaternion<T>& a) {
	return a.sqrmag();
}

template<typename T>
const Vector3<T> ln(const Quaternion<T> &qt) {
	T sc = qt._v.sqrmag();
	T theta = atan2(sc, qt._w);

	if(sc > EPSILON)
		sc = theta / sc;
	else  sc = 1.0 ;

	return sc * qt._v;
}

template<typename T>
const Quaternion<T> exp(const Vector3<T> &v) {
	T theta = mag(v), sc;

	if(theta < EPSILON) sc = 1.0f;
	else sc = sin(theta) / theta;

	return Quaternion<T>(cos(theta), sc * v);
}

template<typename T>
const Quaternion<T> pow(const Vector3<T> &v, T c) {
	return exp(v * c);
}

template<typename T>
const Quaternion<T> slerp(const Quaternion<T>& a, const Quaternion<T>& b, T c) {	

	if( IsAlmostZero(c) )				return a;
	else if( IsAlmostSame(c, (T)1.0f) )	return b;

	// Quaternion's four elements are always kept normalized.
	// So the inner_product directly means the cosine of their angle.
	Quaternion<T> _b = b;

	// To interpolate with closer one.
	// Whichever quaternion get closer when needed, it will produce the same axis-angle, though not the same quaternion but the quaternion in the opposite side.
	_b.setCloser( a );

	T _dot = dot(a, _b);

	// check if they're same.
	// blind copy will be enough.
	if( IsAlmostSame( _dot, (T)1.0f ) )
		return a;   

	T q = acos( _dot );  
	T sinq = sin(q), sin_1_c_q = sin((1.0f-c)*q), sin_c_q = sin(c*q);
	T wa = sin_1_c_q/sinq, wb = sin_c_q/sinq;

	return Quaternion<T>(wa*a._w + wb*b._w, wa*a._v + wb*b._v);
}

template<typename T>
const T distance(const Quaternion<T>& a, const Quaternion<T>& b) {
	return min( (ln( inverse(a) * b )).mag(), (ln( inverse(a) * -b )).mag() );
}

template<typename T>
const Vector3<T> difference(const Quaternion<T>& a, const Quaternion<T>& b) {
	return ln( inverse(b) * a );
}

};	// namespace arithmetic

#endif