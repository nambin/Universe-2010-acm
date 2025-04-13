#ifndef _ARITHMETIC_EULERANGLE_H_
#define _ARITHMETIC_EULERANGLE_H_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "../Vector/Vector3.h"
#include "EulerAngle_Param.h"
#include "Quaternion.h"
#include "RotationMatrix.h"

#pragma warning(disable:4996)

namespace arithmetic {

#ifdef _DEBUG
#define __ASSERT(c) assert(c)
#else
#define __ASSERT(c) {}
#endif

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// NOW, This is only for EA_XYZ_G !!!!!!!!!! (Maybe)
// Maybe, It works for 3 DOF Euler Angle. But Not Repeated..
// Almost I'm not sure... -.-; But It will work for EA_XYZ EulerAngle
//
// INPUT : Degree !!! (Not Radian)
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
// a: axis, p: parity, r: repeat, f: frame
// i, j, and k will be X, Y, and Z respectively
#define GetEulerParam(mode, i, j, k, p, r, f) {				\
	unsigned m = mode;										\
	f=m&1; m>>=1; r=m&1; m>>=1; p=m&1; m>>=1;				\
	i = EulSafe[m&3]; j = EulNext[i+p]; k = EulNext[i+1-p];	\
} 

template<typename T>
class EulerAngle {
	int 			mode;
	unsigned int	DOF;
	Vector3<T> 		ea;

public:
	//////////////////////////////////////////////////////////////////////////
	// Constructor
	EulerAngle() : mode(EA_XYZ_G), DOF(3) {}
	explicit EulerAngle(int m);	
	EulerAngle(int m, const Vector3<T>& v); 
	EulerAngle(int m, T ang1, T ang2, T ang3);
	EulerAngle(int m, T ang1, T ang2);
	EulerAngle(int m, T ang1);
	EulerAngle(int m, const RotationMatrix<T>& M)	{ set(m, M);		}
	EulerAngle(int m, const Quaternion<T>& q)		{ set(m, q);		}
	EulerAngle(int m, const AxisAngle<T>& aa)		{ set(m, aa);		}
	EulerAngle(int m, const Vector3<T>& ax, T ang)	{ set(m, ax, ang);	}

	//////////////////////////////////////////////////////////////////////////
	// Setter
	void set(const EulerAngle<T> &ea)				{ (*this) = ea;						}
	void set(int m, const RotationMatrix<T>& M);
	void set(int m, const Quaternion<T>& q)			{ set(m, RotationMatrix(q));		}
	void set(int m, const AxisAngle<T>& aa)			{ set(m, RotationMatrix(aa));		}
	void set(int m, const Vector3<T>& ax, T ang)	{ set(m, RotationMatrix(ax, ang));	}

	void setMode(int m);

	//////////////////////////////////////////////////////////////////////////
	// Index operators
	inline T&		operator[](unsigned int i)      ;		
	inline const T  operator[](unsigned int i) const;		

	inline T&		operator()(unsigned int i)      ;		
	inline const T  operator()(unsigned int i) const;	

	const int getMode() const;
	const unsigned int getDOF() const	{ return DOF;	}
	const Vector3<T>& getVector() const { return ea;	}
	void setVector(const Vector3<T>& v) { ea = v;		}

	//////////////////////////////////////////////////////////////////////////
	// Manipulation
	inline EulerAngle& invert();
	inline const EulerAngle getInverse() const;

	//////////////////////////////////////////////////////////////////////////
	// Friends Operator
	template <typename T> friend const EulerAngle<T> operator*(const T c, const EulerAngle<T>& ea);
	template <typename T> friend const EulerAngle<T> operator*(const EulerAngle<T>& ea, const T c);

	template <typename T> friend const Vector3<T> operator*(const EulerAngle<T>& ea, const Vector3<T> &v);

	//////////////////////////////////////////////////////////////////////////
	// IO
	template <typename T> friend std::ostream& operator<<(std::ostream& out, const EulerAngle<T>& ea);
	template <typename T> friend std::istream& operator>>(std::istream& in, EulerAngle<T>& ea);

	//////////////////////////////////////////////////////////////////////////
	// Getter _ Other Types
	const Quaternion<T> getQuaternion() const;		// This is only for EA_XYZ_G 
	const Matrix3<T> getMatrix3() const;
	const RotationMatrix<T> getRotationMatrix() const { return getMatrix3(); }
	const AxisAngle<T> getAxisAngle() const;
	void getUnitAxisAngle(Vector3<T>& ax, T& ang) const;

	const Matrix3<T> getInverseMatrix3() const;
};

template<typename T>
inline EulerAngle<T>::EulerAngle(int m) { 
	setMode(m);
}

template<typename T>
inline EulerAngle<T>::EulerAngle(int m, const Vector3<T>& v) : mode(m), DOF(3), ea(v) {
	assert(m <= 0x17);
}

template<typename T>
inline EulerAngle<T>::EulerAngle(int m, T ang1, T ang2, T ang3) : mode(m), DOF(3), ea(ang1, ang2, ang3) {
	assert(m <= 0x17);
}

template<typename T>
inline EulerAngle<T>::EulerAngle(int m, T ang1, T ang2) : mode(m - 0x18), DOF(2), ea(ang1, ang2, 0.0f) {
	assert(m >= 0x18 && m <= 0x2f);
}

template<typename T>
inline EulerAngle<T>::EulerAngle(int m, T ang) : mode(m - 0x30), DOF(1), ea(ang, 0.0f, 0.0f) {
	assert(m >= 0x30);
}

//////////////////////////////////////////////////////////////////////////
// Setter
template<typename T>
void EulerAngle<T>::set(int m, const RotationMatrix<T>& M) {
	assert (m <= 0x17);
	mode = m;
	DOF = 3;

	int i, j, k;
	int p, r, f;
	GetEulerParam(mode, i, j, k, p, r, f); 

	if (r) {
		T sintheta=sqrt(M(i,j)*M(i,j)+M(i,k)*M(i,k));
		if(!basis::IsAlmostZero(sintheta)){
			ea[0] = atan2(M(i,j), M(i,k));
			ea[1] = atan2(sintheta, M(i,i));
			ea[2] = atan2(M(j,i), -M(k,i));
		} 
		else {        // gimbal lock case
			ea[0] = atan2(-M(j,k),M(j,j));
			ea[1] = atan2(sintheta, M(i,i));
			ea[2] = 0.0f;
		}
	} else {
		T costheta=sqrt(M(k,j)*M(k,j)+M(k,k)*M(k,k));
		if(!basis::IsAlmostZero(costheta)){
			ea[0] = atan2(M(k,j), M(k,k));
			ea[1] = atan2(-M(k,i), costheta);
			ea[2] = atan2(M(j,i), M(i,i));
		} 
		else {        // gimbal lock case
			ea[0] = atan2(-M(j,k),M(j,j));
			ea[1] = atan2(-M(k,i), costheta);
			ea[2] = 0.0f;
		}
	}
	if (p) ea = -ea;
	if (f) { std::swap(ea[0], ea[2]); }

	// Degree
	ea *= RadToDeg;
}

template<typename T>
void EulerAngle<T>::setMode(int m) { 
	ea[0] = ea[1] = ea[2] = 0.0f;
	if (m <= 0x17)		{ mode = m;		 DOF = 3;	}
	else if (m <= 0x2f) { mode = m-0x18; DOF = 2;	}
	else				{ mode = m-0x30; DOF = 1;	}
}

//////////////////////////////////////////////////////////////////////////
// Index operators
template <typename T> inline T &EulerAngle<T>::operator[](unsigned int i) {
	__ASSERT(i>=0 && i<DOF); 
	return ea[i];
}
template <typename T> inline T &EulerAngle<T>::operator()(unsigned int i) {
	__ASSERT(i>=0 && i<DOF); 
	return ea[i];
}

template <typename T> inline const T EulerAngle<T>::operator[](unsigned int i) const {
	__ASSERT(i>=0 && i<DOF); 
	return ea[i];
}

template <typename T> inline const T EulerAngle<T>::operator()(unsigned int i) const {
	__ASSERT(i>=0 && i<DOF); 
	return ea[i];
}

//////////////////////////////////////////////////////////////////////////
// Getter
template<typename T>
const int EulerAngle<T>::getMode() const { 
	if(DOF==3)		return mode; 
	else if(DOF==2) return mode+0x18;
	else if(DOF==1) return mode+0x30;

	assert(false);	return 0;
}

//////////////////////////////////////////////////////////////////////////
// Manipulation
template<typename T>
inline EulerAngle<T>& EulerAngle<T>::invert() {
	return (*this) = getInverse();
}

template<typename T>
inline const EulerAngle<T> EulerAngle<T>::getInverse() const {
	return EulerAngle<T>(mode, getMatrix3().getInverse());
}

//////////////////////////////////////////////////////////////////////////
// Friends operators
template <typename T> 
const EulerAngle<T> operator*(const T c, const EulerAngle<T>& EA) {
	if(EA.mode <= 0x17) 	 return EulerAngle(EA.mode, c*EA.ea[0], c*EA.ea[1], c*EA.ea[2]);
	else if(EA.mode <= 0x2f) return EulerAngle(EA.mode, c*EA.ea[0], c*EA.ea[1]);
	else					 return EulerAngle(EA.mode, c*EA.ea[0]);
}

template <typename T> 
const EulerAngle<T> operator*(const EulerAngle<T>& EA, const T c) {
	return c * EA;
}    

template <typename T> 
const Vector3<T> operator*(const EulerAngle<T>& ea, const Vector3<T> &v) {
	return RotationMatrix<T>(ea) * v;
}

//////////////////////////////////////////////////////////////////////////
// IO
template <typename T> 
std::ostream& operator<<(std::ostream& out, const EulerAngle<T>& ea) {
	out << "(" << ea.mode << ", " << ea[0] << ", " << ea[1] << ", " << ea[2] << ")";
	return out;
}

template <typename T> 
std::istream& operator>>(std::istream& in, EulerAngle<T>& EA) {
	char buffer;
	int mode;
	in >>buffer>>mode>>buffer>>EA.ea[0]>>buffer>>EA.ea[1]>>buffer>>EA.ea[2]>>buffer;
	EA.mode = mode;
	return in;
}

//////////////////////////////////////////////////////////////////////////
// Getter _ Other Types
template<typename T>
const Quaternion<T> EulerAngle<T>::getQuaternion() const {				// This is only for EA_XYZ_G 	
	T cr, cp, cy, sr, sp, sy, cpcy, spsy;	
	Vector3<T> EA(ea);
	EA *= static_cast<T>(DegToRad);

	// Formerly roll, pitch, yaw	
	cr = cos(0.5f*EA[0]); cp = cos(0.5f*EA[1]);	cy = cos(0.5f*EA[2]);	
	sr = sin(0.5f*EA[0]); sp = sin(0.5f*EA[1]);	sy = sin(0.5f*EA[2]);	
	cpcy = cp * cy;	
	spsy = sp * sy;	
	return Quaternion<T>(	cr * cpcy + sr * spsy, 
							sr * cpcy - cr * spsy,	
							cr * sp * cy + sr * cp * sy,		
							cr * cp * sy - sr * sp * cy );
}

template<typename T>
const Matrix3<T> EulerAngle<T>::getMatrix3() const {
	int i, j, k;
	int p, r, f;
	GetEulerParam(mode, i, j, k, p, r, f);

	Vector3<T> EA(ea);
	EA *= static_cast<T>(DegToRad);

	if (f) { std::swap(EA[0], EA[2]);	}
	if (p) { EA *= (-1.0f);				}

	T c1 = cos(EA[0]);	T c2 = cos(EA[1]);	T c3 = cos(EA[2]);
	T s1 = sin(EA[0]);	T s2 = sin(EA[1]);	T s3 = sin(EA[2]);

	Matrix3<T> M;
	if (r) {
		M(i,j) = c2;		M(i,j) = s2*s1;				M(i,k) = s2*c1;
		M(j,i) = s2*s3;		M(j,j) = -c2*s1*s3+c1*c3;	M(j,k) = -c2*c1*s3;
		M(k,i) = -s2*c3;	M(k,j) = c2*s1*c3+c1*s3;	M(k,k) = c2*c1*c3-s1*s3;
	} 
	else {
		M(i,i) = c3*c2;		M(i,j) = c3*s2*s1-s3*c1;	M(i,k) = c3*s2*c1+s3*s1;
		M(j,i) = s3*c2;		M(j,j) = s3*s2*s1+c3*c1;	M(j,k) = s3*s2*c1-c3*s1;
		M(k,i) = -s2;		M(k,j) = c2*s1;				M(k,k) = c2*c1;
	}

	return M;
}

template<typename T>
const Matrix3<T> EulerAngle<T>::getInverseMatrix3() const {
	EulerAngle<T> eulerAngle(EA_ZYX_G, -Vector3<T>(ea[2], ea[1], ea[0]));
	return eulerAngle.getMatrix3();
}

template<typename T>
const AxisAngle<T> EulerAngle<T>::getAxisAngle() const {
	Vector3<T> ax; T ang;
	getUnitAxisAngle(ax, ang);
	return AxisAngle(ax, ang);
}

template<typename T>
void EulerAngle<T>::getUnitAxisAngle(Vector3<T>& ax, T& ang) const {
	getQuaternion().getUnitAxisAngle(ax, ang);
}

#undef __ASSERT

//////////////////////////////////////////////////////////////////////////
// Other Class Functions
template<typename T>
inline void Quaternion<T>::set(const EulerAngle<T> &ea) {
	(*this) = ea.getQuaternion();
}

template<typename T>
inline void RotationMatrix<T>::set(const EulerAngle<T> &ea) {
	(*this) = ea.getMatrix3();
}

};		// arithmetic

#endif

