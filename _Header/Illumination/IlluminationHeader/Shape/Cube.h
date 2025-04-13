#ifndef _ILLUMINATION_SHAPE_CUBE_H_
#define _ILLUMINATION_SHAPE_CUBE_H_

#include "Shape.h"

namespace illuminationspace {

class ILLUMINATIONDLL Cube : public ShapeSet {
protected : 
	class CubeValue;

public : 
	Cube(const Transform &o2w, const Vector3<real>& v0, const Vector3<real>& v1)
	: ShapeSet(new CubeValue(o2w, v0, v1)) {}

protected : 
	class ILLUMINATIONDLL CubeValue : public ShapeSet::ShapeSetValue {
	public:	
		//////////////////////////////////////////////////////////////////////////
		// Constructor/Destructor
		CubeValue(const Transform &o2w, const Vector3<real>& v0, const Vector3<real>& v1);

		CubeValue * clone() const { return new CubeValue(*this); }

		//////////////////////////////////////////////////////////////////////////
		// Bounding Box
		const BBox objectBound() const;

		//////////////////////////////////////////////////////////////////////////
		// Manipulator		
		void refine(std::vector<Shape> &refined) const;

		//////////////////////////////////////////////////////////////////////////
		// Property
		const bool canIntersect() const;
		const real area() const	{ return _area;	}

		//////////////////////////////////////////////////////////////////////////
		// Sampling
		const Vector3<real> sample(real u1, real u2, Vector3<real> *Ns) const;

		//////////////////////////////////////////////////////////////////////////
		// Intersection Test
		const bool hit(const Ray &ray, real &t_hit, DifferentialGeometry &dg) const;
		const bool shadowHit(const Ray &ray) const;

	public : 
		BBox cube;
	};
};

};		// illuminationspace

#endif
