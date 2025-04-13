#ifndef _ARITHMETIC_AXISANGLE_H_
#define _ARITHMETIC_AXISANGLE_H_

#include <iostream>
#include "../../Principle.h"
#include "../Vector/Vector3.h"
#include "RotationMatrix.h"
#include "EulerAngle.h"

namespace arithmetic {

using arithmetic::Vector3;

template<typename T>
class AxisAngle	{
private:
	Vector3<T> axis;
	T angle;					// Radian

public:
	inline AxisAngle() : axis(1,0,0), angle(0) {}
	inline AxisAngle(const Vector3<T>& ax, const T ang) : axis(unitVector(ax)), angle(ang) {}

	explicit AxisAngle(const Vector3<T>& aa)		{ set(aa);		}
	AxisAngle(T x, T y, T z)						{ set(x,y,z);	}
	explicit AxisAngle(const RotationMatrix<T>& rt)	{ set(rt);		}
	explicit AxisAngle(const Quaternion<T>& qt)		{ set(qt);		}
	explicit AxisAngle(const EulerAngle<T>& ea)		{ set(ea);		}

	//////////////////////////////////////////////////////////////////////////
	// Setter
	void set(const Vector3<T>& ax, const T ang)	{ setAxis(ax); setAngle(ang);	}
	void setAxis(const Vector3<T>& a)			{ axis = unitVector(a);			}
	void setAngle(const T a)					{ angle = a;					}

	void set(const Vector3<T>& aa)				{ angle = aa.l2Norm(); axis = aa / angle;	}
	void set(T x, T y, T z)						{ set(Vector3<T>(x,y,z));					}

	void set(const RotationMatrix<T>&);
	void set(const Quaternion<T>&);
	void set(const EulerAngle<T>&);	

	//////////////////////////////////////////////////////////////////////////
	// Getter
	const Vector3<T>& getAxis() const		{ return axis;			}
	const T getAngle() const				{ return angle;			}
	const Vector3<T> getAxisAngle() const	{ return angle*axis;	}

	const Quaternion<T> getQuaternion() const;
	const Matrix3<T> getMatrix3() const;

	//////////////////////////////////////////////////////////////////////////
	// Manipulation
	inline AxisAngle& invert()					{ return (*this) = getInverse();	}
	inline const AxisAngle& getInverse() const	{ return AxisAngle<T>(axis, -angle);}

	//////////////////////////////////////////////////////////////////////////
	// Operation
	template<typename T> friend const AxisAngle<T> operator*(const AxisAngle<T>& aa1, const AxisAngle<T>& aa2);

	//////////////////////////////////////////////////////////////////////////
	// IO
	template<typename T> friend std::ostream& operator<<(std::ostream&, const AxisAngle<T>&);
	template<typename T> friend std::istream& operator>>(std::istream&, AxisAngle<T>&);
};

//////////////////////////////////////////////////////////////////////////
// Setter
template<typename T>
void AxisAngle<T>::set(const RotationMatrix<T>& R) {
	T tr = R.trace();

	if(basis::IsAlmostSame(tr, static_cast<T>(3.0f))) {
		angle = 0.0f;
		axis = Vector3<T>(1,0,0);
	}
	else if(basis::IsAlmostSame(tr,static_cast<T>(-1.0f))) {
		angle = static_cast<T>(PI);
		switch (R.indexOfMaxDiagonalComponent()) {
			case 0: axis[0] = sqrt((1.0f + R(0,0)) * 0.5f);
					axis[1] = (R(0,1) + R(1,0)) / (4.0f * axis[0]);
					axis[2] = (R(0,2) + R(2,0)) / (4.0f * axis[0]);
					break;
			case 1: axis[1] = sqrt((1.0f + R(1,1)) * 0.5f);
					axis[2] = (R(1,2) + R(2,1)) / (4.0f * axis[1]);
					axis[0] = (R(1,0) + R(0,1)) / (4.0f * axis[1]);
					break;
			case 2: axis[2] = sqrt((1.0f + R(2,2)) * 0.5f);
					axis[0] = (R(2,0) + R(0,2)) / (4.0f * axis[2]);
					axis[1] = (R(2,1) + R(1,2)) / (4.0f * axis[2]);
					break;
		}
	}
	else {
		T sine = 0.5f * sqrt((R(2,1)-R(1,2))*(R(2,1)-R(1,2)) +
							 (R(0,2)-R(2,0))*(R(0,2)-R(2,0)) +
							 (R(1,0)-R(0,1))*(R(1,0)-R(0,1)) );
		T cosine = (tr - 1.0f) * 0.5f;

		angle = atan2(sine, cosine);
		axis[0] = (R(2,1) - R(1,2)) / (2.0f * sine);
		axis[1] = (R(0,2) - R(2,0)) / (2.0f * sine);
		axis[2] = (R(1,0) - R(0,1)) / (2.0f * sine);
	}
}

template<typename T>
void AxisAngle<T>::set( const Quaternion<T>& qt )	{
	qt.getUnitAxisAngle(axis, angle);
}

template<typename T>
void AxisAngle<T>::set(const EulerAngle<T>& ea) {
	Quaternion qt = ea.getQuaternion();
	qt.getUnitAxisAngle(axis, angle);
}

//////////////////////////////////////////////////////////////////////////
// Getter
template<typename T>
const Quaternion<T> AxisAngle<T>::getQuaternion() const {
	return Quaternion<T>(getAxis(), getAngle());
}

template<typename T>
const Matrix3<T> AxisAngle<T>::getMatrix3() const {
	Vector3<T> ax = getAxis();
	T ang = getAngle();

	Matrix3<T> A = ax.getTensorProduct();	// = ax * trans(ax);
	Matrix3<T> B = (1.0f - cos(ang)) * A;
	Matrix3<T> C = cos(ang) * Matrix3<T>();
	Matrix3<T> D = ax.getDualMatrix();
	Matrix3<T> E = sin(ang) * D;

	return B + C + E;
}

//////////////////////////////////////////////////////////////////////////
// Operation
template<typename T>
inline const AxisAngle<T> operator*(const AxisAngle<T>& aa1, const AxisAngle<T>& aa2) {
	assert(	basis::IsAlmostSame(l2Norm(aa1.getAxis()),(T)1.0f) && 
			basis::IsAlmostSame(l2Norm(aa2.getAxis()),(T)1.0f));

	AxisAngle<T> ret;

	Vector3<T> u1(aa1.getAxis()), u2(aa2.getAxis());
	T theta1 = aa1.getAngle(), theta2 = aa2.getAngle(); 

	T s1 = sin(theta1 * 0.5f);
	T c1 = cos(theta1 * 0.5f);
	T s2 = sin(theta2 * 0.5f);
	T c2 = cos(theta2 * 0.5f);
	Vector3<T> q = s1*c2*u1 + c1*s2*u2 - s1*s2*(cross(u1,u2));
	T nrm = l2Norm(q);

	if(basis::IsAlmostZero(nrm))	ret.setAxis(u1);
	else							ret.setAxis(q * (1.0f/nrm));

	T cosine = c1 * c2 - s1 * s2 * dot(u1,u2);
	T half_theta = atan2(nrm, cosine);

	ret.setAngle((2.0f * half_theta));

	return ret;
}

//////////////////////////////////////////////////////////////////////////
// IO
template<typename T>
std::ostream &operator<<( std::ostream &os, const AxisAngle<T> &val ) {
	os << "(" << val.angle << "," << val.axis[0] << "," << val.axis[1] << "," << val.axis[2] << ")" ;
	return os;
}

template<typename T>
std::istream &operator>>( std::istream &in, AxisAngle<T> &val ) {
	char buffer;

	in >> buffer;	// "("
	in >> val.angle;
	in >> buffer;	// ","
	in >> val.axis[0];	
	in >> buffer;	// ","
	in >> val.axis[1];	
	in >> buffer;	// ","
	in >> val.axis[2];	
	in >> buffer;	// ")"

	val.axis.normalize();
	return in;
}

//////////////////////////////////////////////////////////////////////////
// Other Class Functions
template<typename T> 
void RotationMatrix<T>::set(const AxisAngle<T> &aa) {
	(*this) = aa.getMatrix3();
}

template<typename T> 
inline void Quaternion<T>::set(const AxisAngle<T> &aa) {
	(*this) = aa.getQuaternion();
}

};		// arithmetic

#endif
