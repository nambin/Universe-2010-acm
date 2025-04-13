#ifndef _RIGIDSPACE_RIGIDBODY_H_
#define _RIGIDSPACE_RIGIDBODY_H_

#include "../../../Principle/Arithmetic/Orientation/Orientation.h"
#include "../../../Principle/Geometry/GeometricPrimitive/GeometricObject.h"
#include "../../Nature/DynamicObject.h"

namespace rigidspace {

using nature::DynamicObject;
using geometryspace::GeometricObject;

template<typename VEC>
class RigidBody : public DynamicObject<VEC> {
protected : 
	template<typename VEC> class RigidBodyValue;

public : 	
	RigidBody(RigidBodyValue<VEC> *ptr) : DynamicObject(ptr) {}

	void setFrictionCoefficient(T m) { getptr()->setFrictionCoefficient(m); }

private : 
	RigidBodyValue<VEC>* getptr() { 
		return static_cast<RigidBodyValue<VEC>*>(const_cast<DynamicObjectValue<VEC>*>(value.getptr())); 
	}

protected : 
	template<typename VEC>
	class RigidBodyValue : public DynamicObjectValue<VEC> {
	public : 
		typedef typename VEC::value_type T;

		RigidBodyValue(const GeometricObject<VEC> obj) : object(obj), mu(0) {}

		void		setFrictionCoefficient(T m)					 { mu = m;		}
		const T		getFrictionCoefficient(const VEC& pos) const { return mu;	}

		//////////////////////////////////////////////////////////////////////////
		// glDraw
		void		glDraw() const { object.glDraw(); }

	public : 		
		const GeometricObject<VEC>	object;

		VEC							com;
		Quaternion<T>				qt;	

		// Friction
		T	mu;
	};
};

};		// rigidspace

#endif
