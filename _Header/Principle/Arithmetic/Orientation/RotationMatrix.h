#ifndef _ARITHMETIC_ROTATIONMATRIX_H_
#define _ARITHMETIC_ROTATIONMATRIX_H_

#include "../../Principle.h"
#include "Quaternion.h"

namespace arithmetic {

template<typename T> class AxisAngle;

#ifdef _DEBUG
#define __ASSERT(c) assert(c)
#else
#define __ASSERT(c) {}
#endif

using arithmetic::Vector3;
using arithmetic::Matrix3;

template<typename T>
class RotationMatrix : public Matrix3<T> {
public:	
	RotationMatrix();
	RotationMatrix(const Matrix3<T>& mat);
	RotationMatrix(const Matrix4<T>& mat);
	RotationMatrix(const Vector3<T>& axis, T angle)		{ set(axis, angle);		}
	explicit RotationMatrix(const AxisAngle<T> &aa)		{ set(aa);				}
	explicit RotationMatrix(const Quaternion<T> &qt)	{ set(qt);				}
	explicit RotationMatrix(const EulerAngle<T> &ea)	{ set(ea);				}

	void set(const Matrix3<T>& mat)					{ Matrix3<T>::set(mat);	}
	void set(const Vector3<T>& axis, T angle);
	void set(const AxisAngle<T> &aa);				// Implemented in AxisAngle.h
	void set(const Quaternion<T> &qt)				{ set(qt.getMatrix3());	}
	void set(const EulerAngle<T> &ea);				// Implemented in EulerAngle.h

	// methods	
	inline RotationMatrix &invert();
	const RotationMatrix getInverse() const;
};

template<typename T>
inline RotationMatrix<T>::RotationMatrix() : Matrix3<T>() {}

template<typename T>
inline RotationMatrix<T>::RotationMatrix(const Matrix3<T> &mat) : Matrix3<T>(mat) {
	__ASSERT(isOrthoNormal());
}

template<typename T>
inline RotationMatrix<T>::RotationMatrix(const Matrix4<T> &mat) {
	for(unsigned int i=0;i<3;++i) {
		for(unsigned int j=0;j<3;++j) {
			(*this)(i,j) = mat(i,j);
		}
	}		
	__ASSERT(isOrthoNormal());
}

template<typename T>
void RotationMatrix<T>::set(const Vector3<T>& axis, T angle) {
	Vector3<T> ax = unitVector(axis);
	T ang = angle;

	Matrix3<T> A = ax.getTensorProduct();	// = ax * trans(ax);
	Matrix3<T> B = (1.0f - cos(ang)) * A;
	Matrix3<T> C = cos(ang) * Matrix3<T>();
	Matrix3<T> D = ax.getDualMatrix();
	Matrix3<T> E = sin(ang) * D;

	(*this) = B + C + E;
}

template<typename T>
inline RotationMatrix<T>& RotationMatrix<T>::invert() {
	__ASSERT(isOrthoNormal());

	std::swap((*this)(1,0), (*this)(0,1));
	std::swap((*this)(2,0), (*this)(0,2));
	std::swap((*this)(1,2), (*this)(2,1));

	return (*this);
}

template<typename T>
inline const RotationMatrix<T> RotationMatrix<T>::getInverse() const {
	RotationMatrix<T> R(*this);
	return R.invert();
}

#undef __ASSERT

};		// arithmetic

#endif

