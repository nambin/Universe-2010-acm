#ifndef _RIGIDSPACE_DYNAMIC_RIGIDBODY_H_
#define _RIGIDSPACE_DYNAMIC_RIGIDBODY_H_

#include "RigidBody.h"

namespace rigidspace {

template<typename VEC>
class DynamicRigidBody : public RigidBody<VEC> {
protected : 
	template<typename VEC> class DynamicRigidBodyValue;

public : 
	DynamicRigidBody(const GeometricObject<VEC> obj) : RigidBody(new DynamicRigidBodyValue<VEC>(obj)) {}

protected : 
	template<typename VEC>
	class DynamicRigidBodyValue : public RigidBodyValue<VEC> {
	public : 
		DynamicRigidBodyValue(const GeometricObject<VEC> obj) : RigidBodyValue(obj) {}

		DynamicRigidBodyValue *clone() const { return new DynamicRigidBodyValue(*this); }

		//Vector3<T>	getVel(const Vector3<T>& pos) const;

		////////////////////////////////////////////////////////////////////////////
		//// Distance Field Library
		//bool		isInside(const Vector3<T>& pos) const;

		//T			getPhi(const Vector3<T>& pos) const;
		//Vector3<T>	getNormal(const Vector3<T>& pos) const;

	public : 
		// Dynamic States
		VEC		vel;
		VEC		w;

		// Mass, Inertia
	};
};

};		// rigidspace

#endif
