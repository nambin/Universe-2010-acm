#ifndef _RIGIDSPACE_STATIC_RIGIDBODY_H_
#define _RIGIDSPACE_STATIC_RIGIDBODY_H_

#include "RigidBody.h"

namespace rigidspace {

template<typename VEC>
class StaticRigidBody : public RigidBody<VEC> {
protected : 
	template<typename VEC> class StaticRigidBodyValue;

public : 
	StaticRigidBody(const GeometricObject<VEC> obj) : RigidBody(new StaticRigidBodyValue<VEC>(obj)) {}

protected : 
	template<typename VEC>
	class StaticRigidBodyValue : public RigidBodyValue<VEC> {
	public : 
		StaticRigidBodyValue(const GeometricObject<VEC> obj) : RigidBodyValue(obj) {}

		StaticRigidBodyValue* clone() const { return new StaticRigidBodyValue(*this); }

		const VEC	getVel(const VEC& pos) const		{ return VEC(0);				}

		//////////////////////////////////////////////////////////////////////////
		// Distance Field Library
		const bool	isInside(const VEC& pos) const		{ return object.isInside(pos);	}

		const T		getPhi(const VEC& pos) const		{ return object.getPhi(pos);	}
		const VEC	getNormal(const VEC& pos) const		{ return object.getNormal(pos); }
	};
};

template<typename T>
class StaticRigidBody2D : public StaticRigidBody<Vector2<T>> {
public : 
	StaticRigidBody2D(const GeometricObject<Vector2<T>> obj) : StaticRigidBody<Vector2<T>>(obj) {}
};

template<typename T>
class StaticRigidBody3D : public StaticRigidBody<Vector3<T>> {
public : 
	StaticRigidBody3D(const GeometricObject<Vector3<T>> obj) : StaticRigidBody<Vector3<T>>(obj) {}
};

};		// rigidspace

#endif
