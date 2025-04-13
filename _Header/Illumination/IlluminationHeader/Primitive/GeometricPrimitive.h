#ifndef _ILLUMINATION_GEOMETRIC_PRIMITIVE_H_
#define _ILLUMINATION_GEOMETRIC_PRIMITIVE_H_

#include "Primitive.h"

namespace illuminationspace {

class ILLUMINATIONDLL GeometricPrimitive : public Primitive {
protected : 
	class GeometricPrimitiveValue;

public : 
	GeometricPrimitive(const Shape& s, const Material &m, AreaLight a = AreaLight())
	: Primitive(new GeometricPrimitiveValue(s,m,a)) {}

protected : 
	class ILLUMINATIONDLL GeometricPrimitiveValue : public Primitive::PrimitiveValue {
	public : 
		GeometricPrimitiveValue(const Shape &s, const Material &m, AreaLight a = AreaLight())
		: shape(s), material(m), areaLight(a) {}

		GeometricPrimitiveValue * clone() const { return new GeometricPrimitiveValue(*this); }

		const BBox worldBound() const	{ return shape.worldBound();	}
		const bool canIntersect() const { return shape.canIntersect();	}
		void refine(std::vector<Primitive> &refined) const;

		const AreaLight& getAreaLight() const				{ return areaLight; }
		BSDF *getBSDF(const DifferentialGeometry &dg) const	{ return material.getBSDF(dg); }
		const Spectrum extinction(const DifferentialGeometry& dg, const Vector3<real>& wi) const { return material.extinction(dg, wi); }

		const bool hit(const Ray& r, Intersection &intersection) const;
		const bool shadowHit(const Ray& r) const;

	public : 
		Shape		shape;
		Material	material;
		AreaLight	areaLight;
	};	
};

};		// illuminationspace

#endif