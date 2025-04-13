#ifndef _ILLUMINATION_SPHERE_H_
#define _ILLUMINATION_SPHERE_H_

#include "Shape.h"

namespace illuminationspace {

class ILLUMINATIONDLL Sphere : public Shape {
public :
	Sphere(const Transform &o2w, const Vector3<real>& _center, real _radius) 
	: Shape(new SphereValue(o2w, _center, _radius)) {}

protected : 
	class ILLUMINATIONDLL SphereValue : public Shape::ShapeValue {
	public : 
		SphereValue(const Transform &o2w, const Vector3<real>& _center, real _radius) 
		: Shape::ShapeValue(o2w), center(_center), radius(_radius) {}

		SphereValue *clone() const { return new SphereValue(*this); }

		const BBox objectBound() const;

		const real area() const { return 4.0f * PI_i * radius * radius; }

		const bool hit(const Ray& ray, real &tHit, DifferentialGeometry& geometry) const;
		const bool shadowHit(const Ray& ray) const;	

	public :
		Vector3<real> center;
		real radius;
	};
};

};		// illuminationspace

#endif

