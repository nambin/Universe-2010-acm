#ifndef _ILLUMINATION_ISOSURFACE_H_
#define _ILLUMINATION_ISOSURFACE_H_

#include "Shape.h"

namespace illuminationspace {

class ILLUMINATIONDLL IsoSurface : public Shape {
protected : 
	class IsoSurfaceValue;

public:
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	IsoSurface(IsoSurfaceValue *ptr) : Shape(ptr) {}
	
protected : 
	class ILLUMINATIONDLL IsoSurfaceValue : public Shape::ShapeValue {
	public : 
		IsoSurfaceValue() : Shape::ShapeValue() {}
		IsoSurfaceValue(const Transform& o2w) : Shape::ShapeValue(o2w) {}

		const real secant_iteration(const real v[8], const Vector3<real>& w0, const Vector3<real>& w_icr, const real tol, const bool need_exact_time) const;

	protected : 
		const real secant_iteration(const real v[4], const real tol, const bool need_exact_time) const;
	};
};

};		// illuminationspace

#endif