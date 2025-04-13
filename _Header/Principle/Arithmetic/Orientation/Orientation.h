#ifndef _ARITHMETIC_ORIENTATION_H_
#define _ARITHMETIC_ORIENTATION_H_

#include "AxisAngle.h"
#include "EulerAngle.h"
#include "Quaternion.h"
#include "RotationMatrix.h"

namespace arithmetic {

#define AXISANGLE 		(1 << 0)
#define ROTATIONMATRIX	(1 << 1)
#define EULERANGLE		(1 << 2)
#define QUATERNION		(1 << 3)

template<typename T>
class Orientation {
private:
	char				state;
	RotationMatrix<T> 	rm;
	EulerAngle<T>		ea;
	Quaternion<T>		qt;
	AxisAngle<T>		aa;

public:
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Orientation()								{ identity();	}
	Orientation(const RotationMatrix<T>& src)	{ set(src);		}
	Orientation(const EulerAngle<T>& src)		{ set(src);		}
	Orientation(const Quaternion<T>& src)		{ set(src);		}
	Orientation(const AxisAngle<T>& src)		{ set(src);		}

	//////////////////////////////////////////////////////////////////////////
	// Setter
	void identity()								{ set(RotationMatrix<T>());				}
	void set(const RotationMatrix<T>& src)		{ state = ROTATIONMATRIX;	rm = src;	}
	void set(const EulerAngle<T>& src)			{ state = EULERANGLE;		ea = src;	}
	void set(const Quaternion<T>& src)			{ state = QUATERNION;		qt = src;	}
	void set(const AxisAngle<T>& src)			{ state = AXISANGLE;		aa = src;	}

	void setEulerAngle(const Vector3<T>& v)		{ state = EULERANGLE; ea.setVector(v);	}
	void setAxisAngle(const Vector3<T>& v)		{ state = AXISANGLE;  aa.set(v);		}
	void setAxisAngle(T x, T y, T z)			{ setAxisAngle(Vector3<T>(x,y,z));		}

	//////////////////////////////////////////////////////////////////////////
	// Getter
	const RotationMatrix<T>& getRotationMatrix() const;
	const EulerAngle<T>& getEulerAngle() const; 
	const Quaternion<T>& getQuaternion() const;
	const AxisAngle<T>& getAxisAngle() const;

	//////////////////////////////////////////////////////////////////////////
	// Operators
	Orientation& operator*=(const Orientation&);
	template<typename T> friend const Orientation<T> operator*(const Orientation<T>& a, const Orientation<T> &b);
	
	//////////////////////////////////////////////////////////////////////////
	// Manipulation
	Orientation& invert();
	const Orientation getInverse() const;
};

template<typename T>
const RotationMatrix<T>& Orientation<T>::getRotationMatrix() const {
	if(state & ROTATIONMATRIX) return rm;

	state |= ROTATIONMATRIX;

	if(state & EULERANGLE)		rm = ea.getMatrix3();
	else if(state & QUATERNION) rm = qt.getMatrix3();
	else if(state & AXISANGLE)	rm = aa.getMatrix3();
	else 
		assert(0); // Something wrong in orientation.

	return rm;
}

template<typename T>
const EulerAngle<T>& Orientation<T>::getEulerAngle() const {
	if(state & EULERANGLE) return ea;

	state |= EULERANGLE;

	if(state & ROTATIONMATRIX)	ea = EulerAngle<T>(ea.getMode(), rm);
	else if(state & QUATERNION) ea = EulerAngle<T>(ea.getMode(), qt);
	else if(state & AXISANGLE)	ea = EulerAngle<T>(ea.getMode(), aa);
	else 
		assert(0); // Something wrong in orientation.

	return ea;
}

template<typename T>
const Quaternion<T>& Orientation<T>::getQuaternion() const {
	if(state & QUATERNION) return qt;

	state |= QUATERNION;

	if(state & ROTATIONMATRIX)	qt = Quaternion<T>(rm);
	else if(state & EULERANGLE) qt = ea.getQuaternion();
	else if(state & AXISANGLE)	qt = aa.getQuaternion();
	else 
		assert(0); // Something wrong in orientation.

	return qt;
}

template<typename T>
const AxisAngle<T>& Orientation<T>::getAxisAngle() const {
	if(state & AXISANGLE) return aa;

	state |= AXISANGLE;

	if(state & ROTATIONMATRIX)	aa = AxisAngle<T>(rm);
	else if(state & EULERANGLE) aa = AxisAngle<T>(ea);
	else if(state & QUATERNION) aa = AxisAngle<T>(qt);
	else 
		assert(0); // Something wrong in orientation.

	return aa;
}

template<typename T>
inline Orientation<T>& Orientation<T>::operator*=(const Orientation<T> &ori) {
	return (*this) = (*this) * ori;
}

template<typename T>
inline const Orientation<T> operator*(const Orientation<T>& o1, const Orientation<T>& o2) {
	return Orientation<T>(o1.getRotationMatrix() * o2.getRotationMatrix());
}

template<typename T>
Orientation<T>& Orientation<T>::invert() {
	return (*this) = getInverse();
}

template<typename T>
const Orientation<T> Orientation<T>::getInverse() const {
	if(state & QUATERNION)		return Orientation(qt.getInverse());
	if(state & AXISANGLE)		return Orientation(aa.getInverse());
	if(state & ROTATIONMATRIX)	return Orientation(rm.getInverse());
	if(state & EULERANGLE)		return Orientation(ea.getInverse());
	assert(0); return (*this);
}

};		// arithmetic

#endif