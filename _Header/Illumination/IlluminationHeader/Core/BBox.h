#ifndef _ILLUMINATION_BBOX_H_
#define _ILLUMINATION_BBOX_H_

#include "../../IlluminationSpace.h"
#include "Ray.h"
#include "../../../Principle/Geometry/GeometricPrimitive/AABB.h"

namespace illuminationspace {

class ILLUMINATIONDLL BBox : public geometryspace::AABB<real> {
public : 
	static inline const real Large()	{ return LARGE_i; }

public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructors
	inline BBox() : geometryspace::AABB<real>() {}
	inline BBox(const geometryspace::AABB<real>& aabb) : geometryspace::AABB<real>(aabb) {}
	inline BBox(const Vector3<real>& pt) : geometryspace::AABB<real>(pt) {}
	inline BBox(const Vector3<real>& v1, const Vector3<real>& v2) : geometryspace::AABB<real>(v1, v2) {}
	
	//////////////////////////////////////////////////////////////////////////
	// Checker	
	const bool rayIntersect(const Ray& r) const;
	const bool intersectP(const Ray& r, real &t0, real &t1) const;
};

};		// illuminationspace

#endif
