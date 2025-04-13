#include <stdafx.h>
#include <Illumination/IlluminationHeader/Shape/Triangle.h>
#include <Illumination/IlluminationHeader/Core/IlluminationLibrary.h>

using namespace illuminationspace;

const BBox Triangle::TriangleValue::objectBound() const	{
	const Vector3<real> &p0 = meshptr->posVerts(v[0]);
	const Vector3<real> &p1 = meshptr->posVerts(v[1]);
	const Vector3<real> &p2 = meshptr->posVerts(v[2]);

	BBox bbox(WorldToObject.transformPoint(p0));
	bbox.surround(WorldToObject.transformPoint(p1));
	bbox.surround(WorldToObject.transformPoint(p2));

	return bbox;
}

const BBox Triangle::TriangleValue::worldBound() const {

	const Vector3<real> &p0 = meshptr->posVerts(v[0]);
	const Vector3<real> &p1 = meshptr->posVerts(v[1]);
	const Vector3<real> &p2 = meshptr->posVerts(v[2]);

	BBox bbox(p0);
	bbox.surround(p1);
	bbox.surround(p2);

	return bbox;
}

const real Triangle::TriangleValue::area() const {
	const Vector3<real> &p0 = meshptr->posVerts(v[0]);
	const Vector3<real> &p1 = meshptr->posVerts(v[1]);
	const Vector3<real> &p2 = meshptr->posVerts(v[2]);

	return 0.5f * mag(cross(p1-p0, p2-p0));
}	

const Vector3<real> Triangle::TriangleValue::sample(real u1, real u2, Vector3<real> *N) const	{
	real b0, b1;
	IlluminationLibrary::UniformSampleTriangle(u1, u2, &b0, &b1);

	const Vector3<real> &p0 = meshptr->posVerts(v[0]);
	const Vector3<real> &p1 = meshptr->posVerts(v[1]);
	const Vector3<real> &p2 = meshptr->posVerts(v[2]);

	Vector3<real> p = b0 * p0 + b1 * p1 + (1.0f - b0 - b1) * p2;
	*N = uvw->w();

	return p;
}

const bool Triangle::TriangleValue::hit(const Ray& r, real &tHit, DifferentialGeometry &geometry) const {

	const Vector3<real> &p1 = meshptr->posVerts(v[0]);
	const Vector3<real> &p2 = meshptr->posVerts(v[1]);
	const Vector3<real> &p3 = meshptr->posVerts(v[2]);

	Vector3<real> e1 = p2 - p1;
	Vector3<real> e2 = p3 - p1;
	Vector3<real> s1 = arithmetic::cross(r.direction(), e2);

	real divisor = arithmetic::dot(s1, e1);
	if(divisor == 0.0f)
		return false;

	real invDivisor = 1.0f / divisor;

	Vector3<real> d = r.origin() - p1;
	real b1 = arithmetic::dot(d, s1) * invDivisor;
	if(b1 < 0.0f || b1 > 1.0f) return false;

	Vector3<real> s2 = arithmetic::cross(d, e1);
	real b2 = arithmetic::dot(r.direction(), s2) * invDivisor;

	if(b2 < 0.0f || b1 + b2 > 1.0f) return false;

	real t = arithmetic::dot(e2, s2) * invDivisor;
	if(t < r.tMin() || t > r.tMax()) return false;

	tHit = t;
	geometry.hit_p = r(tHit);
	geometry.uvw = (*uvw);

	if(meshptr->normals.size() == 0)
		geometry.shadingNormal = uvw->w();
	else {
		Vector3<real> _normal = (1.0f-b1-b2)*meshptr->normals(n[0]) + b1*meshptr->normals(n[1]) + b2*meshptr->normals(n[2]);
		geometry.shadingNormal = unitVector(_normal);
	}

	if(meshptr->uvs.size() == 0)		
		geometry.uv = Vector2<real>(0.0f);
	else 
		geometry.uv = (1.0f-b1-b2)*meshptr->uvs(v[0]) + b1*meshptr->uvs(v[1]) + b2*meshptr->uvs(v[2]);

	return true;
}

const bool Triangle::TriangleValue::shadowHit(const Ray& r) const {

	const Vector3<real> &p1 = meshptr->posVerts(v[0]);
	const Vector3<real> &p2 = meshptr->posVerts(v[1]);
	const Vector3<real> &p3 = meshptr->posVerts(v[2]);

	Vector3<real> e1 = p2 - p1;
	Vector3<real> e2 = p3 - p1;
	Vector3<real> s1 = arithmetic::cross(r.direction(), e2);

	real divisor = arithmetic::dot(s1, e1);
	if(divisor == 0.0f)
		return false;

	real invDivisor = 1.0f / divisor;

	Vector3<real> d = r.origin() - p1;
	real b1 = arithmetic::dot(d, s1) * invDivisor;
	if(b1 < 0.0f || b1 > 1.0f) return false;

	Vector3<real> s2 = arithmetic::cross(d, e1);
	real b2 = arithmetic::dot(r.direction(), s2) * invDivisor;

	if(b2 < 0.0f || b1 + b2 > 1.0f) return false;

	real t = arithmetic::dot(e2, s2) * invDivisor;
	if(t < r.tMin() || t > r.tMax()) return false;

	return true;
}
