#ifndef _GEOMETRYSPACE_INTERSECTIONTEST_LIBRARY_H_
#define _GEOMETRYSPACE_INTERSECTIONTEST_LIBRARY_H_

#include "../../Principle.h"
#include "../../../Basis/Basis/BasicLibrary.h"
#include "../../Arithmetic/Vector/Vector3.h"

namespace geometryspace {

using arithmetic::Vector3;

class PRINCIPLEDLL IntersectionTestLibrary {
public : 
	template<typename T>
	static const bool testAABBPlane(const Vector3<T> &_min, const Vector3<T> &_max, const Vector3<T> &p_normal, T p_dist);

	template<typename T>
	static const bool testTriangleAABB(const Vector3<T>& v0, const Vector3<T>& v1, const Vector3<T>& v2, const Vector3<T> &_min, const Vector3<T> &_max);	
};

template<typename T>
const bool IntersectionTestLibrary::testAABBPlane(const Vector3<T> &_min, const Vector3<T> &_max, const Vector3<T> &p_normal, T p_dist)	{

	// These two lines not necessary with a (center, extents) AABB representation    
	const Vector3<T> c = (_min + _max) * static_cast<T>(0.5f);	// Compute AABB center
	const Vector3<T> e = _max - c;								// Compute positive extents

	// Compute the projection interval radius of b onto L(t) = b.c + t * p.n
	const T r = e[0] * abs(p_normal[0]) + e[1] * abs(p_normal[1]) + e[2] * abs(p_normal[2]);

	// Compute distance of box center from plane
	const T s = dot(p_normal, c) - p_dist;

	// Intersection occurs when distance s falls within [-r,+r] interval
	return abs(s) <= r;
}

template<typename T>
const bool IntersectionTestLibrary::testTriangleAABB(const Vector3<T>& _v0, const Vector3<T>& _v1, const Vector3<T>& _v2, const Vector3<T> &_min, const Vector3<T> &_max)
{
	T p0, p1, p2, r;

	const Vector3<T> c = (_min + _max) * static_cast<T>(0.5f);	// Compute AABB center

	const T e0 = (_max[0] - _min[0]) * static_cast<T>(0.5f);
	const T e1 = (_max[1] - _min[1]) * static_cast<T>(0.5f);
	const T e2 = (_max[2] - _min[2]) * static_cast<T>(0.5f);

	const Vector3<T> v0 = _v0 - c;
	const Vector3<T> v1 = _v1 - c;
	const Vector3<T> v2 = _v2 - c;

	const Vector3<T> f0 = v1 - v0;
	const Vector3<T> f1 = v2 - v1;
	const Vector3<T> f2 = v0 - v2;

	// Test the three axes corresponding to the face normals of AABB b (category 1)
	if( MAX3(v0[0], v1[0], v2[0]) < -e0 || MIN3(v0[0], v1[0], v2[0]) > e0) return false;
	if( MAX3(v0[1], v1[1], v2[1]) < -e1 || MIN3(v0[1], v1[1], v2[1]) > e1) return false;
	if( MAX3(v0[2], v1[2], v2[2]) < -e2 || MIN3(v0[2], v1[2], v2[2]) > e2) return false;

	// Test axes a00..a22 (category 3)
	// a00
	p0 = v0[2] * v1[1] - v0[1] * v1[2];
	p2 = v2[2] * (v1[1] - v0[1]) - v2[1] * (v1[2] - v0[2]);
	r = e1 * abs(f0[2]) + e2 * abs(f0[1]);
	if(MAX(-MAX(p0, p2), MIN(p0, p2)) > r)
		return false;

	// a01
	p0 = - v0[1] * (v2[2] - v1[2]) + v0[2] * (v2[1] - v1[1]);
	p1 = v1[2] * v2[1] - v1[1] * v2[2];
	r = e1 * abs(f1[2]) + e2 * abs(f1[1]);
	if(MAX(-MAX(p0, p1), MIN(p0, p1)) > r)
		return false;

	// a02
	p0 = v0[1] * v2[2] - v0[2] * v2[1];
	p1 = - v1[1] * (v0[2] - v2[2]) + v1[2] * (v0[1] - v2[1]);
	r = e1 * abs(f2[2]) + e2 * abs(f2[1]);
	if(MAX(-MAX(p0, p1), MIN(p0, p1)) > r)
		return false;

	// a10
	p0 = v0[0] * v1[2] - v0[2] * v1[0];
	p2 = v2[0] * (v1[2] - v0[2]) - v2[2] * (v1[0] - v0[0]);
	r = e0 * abs(f0[2]) + e2 * abs(f0[0]);
	if(MAX(-MAX(p0, p2), MIN(p0, p2)) > r)
		return false;

	// a11
	p0 = v0[0] * (v2[2] - v1[2]) - v0[2] * (v2[0] - v1[0]);
	p1 = v1[0] * v2[2] - v1[2] * v2[0];
	r = e0 * abs(f1[2]) + e2 * abs(f1[0]);
	if(MAX(-MAX(p0, p1), MIN(p0, p1)) > r)
		return false;

	// a12
	p0 = v0[2] * v2[0] - v0[0] * v2[2];
	p1 = v1[0] * (v0[2] - v2[2]) - v1[2] * (v0[0] - v2[0]);
	r = e0 * abs(f2[2]) + e2 * abs(f2[0]);
	if(MAX(-MAX(p0, p1), MIN(p0, p1)) > r)
		return false;

	// a20
	p0 = v0[1] * v1[0] - v0[0] * v1[1];
	p2 = - v2[0] * (v1[1] - v0[1]) + v2[1] * (v1[0] - v0[0]);
	r = e0 * abs(f0[1]) + e1 * abs(f0[0]);
	if(MAX(-MAX(p0, p2), MIN(p0, p2)) > r)
		return false;

	// a21
	p0 = - v0[0] * (v2[1] - v1[1]) + v0[1] * (v2[0] - v1[0]);
	p1 = v2[0] * v1[1] - v1[0] * v2[1];
	r = e0 * abs(f1[1]) + e1 * abs(f1[0]);
	if(MAX(-MAX(p0, p1), MIN(p0, p1)) > r)
		return false;

	// a22
	p0 = v0[0] * v2[1] - v0[1] * v2[0];
	p1 = - v1[0] * (v0[1] - v2[1]) + v1[1] * (v0[0] - v2[0]);
	r = e0 * abs(f2[1]) + e1 * abs(f2[0]);
	if(MAX(-MAX(p0, p1), MIN(p0, p1)) > r)
		return false;

	// Test separating axis corresponding to triangle face normal (category 2)	
	const Vector3<T> p_normal = cross(f0, f1); //p_normal.normalize();
	const T p_dist = arithmetic::dot(p_normal, _v0);

	return testAABBPlane(_min, _max, p_normal, p_dist);
}

};	// geometryspace

#endif
