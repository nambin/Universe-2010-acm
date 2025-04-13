#ifndef _FLUIDSPACE_VELOCITY_FIELDS_H_
#define _FLUIDSPACE_VELOCITY_FIELDS_H_

#include "../../../Principle/Arithmetic/Vector/Vector2.h"
#include "../../Nature.h"
#include "../Environment/Environments.h"

namespace fluidspace {

using arithmetic::Vector2;
using arithmetic::Vector3;

template<typename T>
class VelocityField {
public : 
	VelocityField() { time = 0; }
	virtual ~VelocityField() {}

	virtual const unsigned int getDimension() const = 0;

	void setTime(T t) { time = t; }
	void addTime(T t) { time += t; }

	virtual const T getCFL(T dt, T h) const = 0;

protected : 
	T	time;
};

template<typename T>
class VelocityField2D : public VelocityField<T> {
public : 
	VelocityField2D() : VelocityField() { available_maxVel = false; }

	const unsigned int getDimension() const { return 2; }

	const Vector2<T> operator()(const Vector2<T> &pos) const	 { return getVel(pos);						}
	const Vector2<T> getVel(const Vector2<T> &pos) const		 { return Vector2<T>(getU(pos), getV(pos));	}

	const T getCFL(T dt, T h) const {
		if(available_maxVel == false)
			calc_maxVel();

		available_maxVel = true;
		return dt * (abs(maxVel[0]) + abs(maxVel[1])) / h;
	}
	const Vector2<T> get_maxVel() const { if(available_maxVel == false) calc_maxVel(); return maxVel; }

	virtual const T getU(const Vector2<T> &pos) const = 0;
	virtual const T getV(const Vector2<T> &pos) const = 0;

	mutable bool		available_maxVel;
	mutable Vector2<T>	maxVel;

protected : 
	virtual void calc_maxVel() const = 0;
};

template<typename T>
class VelocityField3D : public VelocityField<T> {
public : 
	VelocityField3D() : VelocityField() { available_maxVel = false; }

	virtual VelocityField3D<T>* clone() const = 0;

	const unsigned int getDimension() const { return 3; }

	const Vector3<T> operator()(const Vector3<T> &pos) const	 { return getVel(pos);									}
	virtual const Vector3<T> getVel(const Vector3<T> &pos) const { return Vector3<T>(getU(pos), getV(pos), getW(pos));	}

	const T getCFL(T dt, T h) const {
		if(available_maxVel == false)
			calc_maxVel();

		available_maxVel = true;
		return dt * (abs(maxVel[0]) + abs(maxVel[1]) + abs(maxVel[2])) / h;
	}

	const Vector3<T> get_maxVel() const { if(available_maxVel == false) calc_maxVel(); return maxVel; }

	virtual const T getU(const Vector3<T> &pos) const = 0;
	virtual const T getV(const Vector3<T> &pos) const = 0;
	virtual const T getW(const Vector3<T> &pos) const = 0;

	mutable bool		available_maxVel;
	mutable Vector3<T>	maxVel;

protected : 
	virtual void calc_maxVel() const = 0;
};

};		// fluidspace

#include "VelocityFields2D.h"
#include "VelocityFields3D.h"
#include "VelocityField_Advector.h"

#endif