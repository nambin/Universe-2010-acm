#ifndef _NATURE_DYNAMIC_OBJECT_H_
#define _NATURE_DYNAMIC_OBJECT_H_

#include "../../../Basis/MemoryManagement/RCPtr.h"
#include "../../../Principle/Arithmetic/Vector/Vector2.h"
#include "../../../Principle/Arithmetic/Vector/Vector3.h"

namespace nature {

template<typename VEC>
class DynamicObject {
protected : 
	template<typename VEC> class DynamicObjectValue;	

public : 
	typedef typename VEC::value_type T;

	DynamicObject() : value(NULL) {}
	DynamicObject(DynamicObjectValue<VEC> *ptr) : value(ptr) {}

	const VEC	getVel(const VEC& pos) const { return value->getVel(pos); }  

	//////////////////////////////////////////////////////////////////////////
	// Distance Field Library
	const bool	isInside(const VEC& pos) const { return value->isInside(pos); }

	const T		getPhi(const VEC& pos) const { return value->getPhi(pos); }
	const VEC	getNormal(const VEC& pos) const { return value->getNormal(pos); }

	const T		getFrictionCoefficient(const VEC& pos) const { return value->getFrictionCoefficient(pos); }

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	void		glDraw() const { value->glDraw(); }

	const unsigned int getDimension() const { return VEC::getDimension(); }

protected : 
	template<typename VEC> 
	class DynamicObjectValue : public basis::RCObject {
	public : 
		typedef typename VEC::value_type T;

		DynamicObjectValue() {}

		virtual DynamicObjectValue * clone() const = 0;	

		virtual const VEC	getVel(const VEC& pos) const = 0; 

		//////////////////////////////////////////////////////////////////////////
		// Distance Field Library
		virtual const bool	isInside(const VEC& pos) const = 0;

		virtual const T		getPhi(const VEC& pos) const = 0;
		virtual const VEC	getNormal(const VEC& pos) const = 0;

		virtual const T		getFrictionCoefficient(const VEC& pos) const = 0;

		//////////////////////////////////////////////////////////////////////////
		// glDraw
		virtual void		glDraw() const = 0;
	};

	//////////////////////////////////////////////////////////////////////////
	// Reference Counted Pointer
	RCPtr<DynamicObjectValue<VEC>> value;

public : 
	const RCPtr<DynamicObjectValue<VEC>>& getDynamicObjectValue() const { return value; }
};

};		// nature

#endif