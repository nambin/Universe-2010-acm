#ifndef _GEOMETRYSPACE_GEOMETRIC_OBJECT_H_
#define _GEOMETRYSPACE_GEOMETRIC_OBJECT_H_

#include <cstdlib>
#ifdef _WIN32
#include "gl/glut.h"
#else
#include <GLUT/glut.h>
#endif

#include "../../../Basis/MemoryManagement/RCPtr.h"
#include "../../Arithmetic/Vector/Vector2.h"
#include "../../Arithmetic/Vector/Vector3.h"
#include "../../VirtualSpace/Material/Material.h"

namespace geometryspace {
	
using virtualspace::Material;
using arithmetic::Vector2;
using arithmetic::Vector3;

template<typename VEC>
class GeometricObject {
protected : 
	template<typename VEC> class GeometricObjectValue;

public : 
	typedef typename VEC::value_type T;

	GeometricObject() : value(NULL) {}
	GeometricObject(GeometricObjectValue<VEC> *ptr) : value(ptr) {}
		
	//////////////////////////////////////////////////////////////////////////
	// Distance Field Library
	const bool	isInside(const VEC& pos) const { return value->isInside(pos); }

	const T		getPhi(const VEC& pos) const { return value->getPhi(pos); }
	const VEC	getNormal(const VEC& pos) const { return value->getNormal(pos); }

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	void		glDraw() const { value->glDraw(); }

protected :
	template<typename VEC>
	class GeometricObjectValue : public basis::RCObject {
	public : 
		GeometricObjectValue() {}
		GeometricObjectValue(const Material& m) : mat(m) {}

		virtual GeometricObjectValue* clone() const = 0;	

		//////////////////////////////////////////////////////////////////////////
		// Distance Field Library
		virtual const bool	isInside(const VEC&) const = 0;

		virtual const T		getPhi(const VEC&) const = 0;
		virtual const VEC	getNormal(const VEC&) const = 0;

		//////////////////////////////////////////////////////////////////////////
		// glDraw
		virtual void		glDraw() const = 0;

	public : 
		Material mat;
	};

	//////////////////////////////////////////////////////////////////////////
	// Reference Counted Pointer
	RCPtr<GeometricObjectValue<VEC>> value;

public : 
	const RCPtr<GeometricObjectValue<VEC>>& getGeometricObjectValue() const { return value; }
};

};		// geometryspace

#endif
