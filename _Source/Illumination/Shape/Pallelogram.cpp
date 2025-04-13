#include <stdafx.h>
#include <Illumination/IlluminationHeader/Shape/Pallelogram.h>

using namespace illuminationspace;

Parallelogram::Parallelogram(const Transform &o2w, const Vector3<real>& _base, const Vector3<real>& _u, const Vector3<real>& _v)
: Shape(new Parallelogram::ParallelogramValue(o2w, _base, _u, _v)) {}

Parallelogram::Parallelogram(const Transform &o2w, const Vector3<real>& _base, const Vector3<real>& _u, const Vector3<real>& _v, 
			  const Vector2<real> & _uv0, const Vector2<real>& _uv1, const Vector2<real>& _uv2)
: Shape(new Parallelogram::ParallelogramValue(o2w, _base, _u, _v, _uv0, _uv1, _uv2)) {}

Parallelogram::ParallelogramValue::ParallelogramValue(const Transform &o2w, const Vector3<real>& _base, const Vector3<real>& _u, const Vector3<real>& _v) 
: ShapeValue(o2w), base(_base), u(_u), v(_v) {

	base = ObjectToWorld.transformPoint(base);
	u = ObjectToWorld.transformVector(u);
	v = ObjectToWorld.transformVector(v);

	norm = unitVector(cross(u,v));
	unorm = unitVector(u);
	vnorm = unitVector(v);

	uvw.initFromUnitVector(norm);

	uv0 = Vector2<real>(0.0f,0.0f);
	uv1 = Vector2<real>(1.0f,0.0f);
	uv2 = Vector2<real>(0.0f,1.0f);

	// find area of Parallelogram
	real width = u.mag();
	Vector3<real> proj = v - ((dot(v,u) / (width*width)) * u);
	real height = proj.mag();

	_area = width * height;
	_pdf = 1.0f / _area;
}

Parallelogram::ParallelogramValue::ParallelogramValue(const Transform &o2w, const Vector3<real> &_base, const Vector3<real> &_u, const Vector3<real> &_v, 
													  const Vector2<real> &_uv0, const Vector2<real> &_uv1, const Vector2<real> &_uv2)
: ShapeValue(o2w), base(_base), u(_u), v(_v), uv0(uv0), uv1(uv1), uv2(uv2) {

	base = ObjectToWorld.transformPoint(base);
	u = ObjectToWorld.transformVector(u);
	v = ObjectToWorld.transformVector(v);

	norm = unitVector(cross(u,v));
	unorm = unitVector(u);
	vnorm = unitVector(v);

	uvw.initFromUnitVector(norm);

	// find area of Parallelogram::ParallelogramValue
	real width = u.mag();
	Vector3<real> proj = v - ((dot(v,u) / (width*width)) * u);
	real height = proj.mag();

	_area = width * height;
	_pdf = 1.0f / _area;
}

const BBox Parallelogram::ParallelogramValue::objectBound() const	{

	Vector3<real> p0 = base;
	Vector3<real> p1 = base + u;
	Vector3<real> p2 = base + v;
	Vector3<real> p3 = base + u + v;

	BBox objectBounds(WorldToObject.transformPoint(p0));		objectBounds.surround(WorldToObject.transformPoint(p1));
	objectBounds.surround(WorldToObject.transformPoint(p2));	objectBounds.surround(WorldToObject.transformPoint(p3));

	return objectBounds;
}

const BBox Parallelogram::ParallelogramValue::worldBound() const	{

	Vector3<real> min;
	Vector3<real> max;

	Vector3<real> p0 = base;
	Vector3<real> p1 = base + u;
	Vector3<real> p2 = base + v;
	Vector3<real> p3 = base + u + v;

	min = minVec(p0, p1);	min = minVec(min, p2);	min = minVec(min, p3);
	max = maxVec(p0, p1);	max = maxVec(max, p2);	max = maxVec(max, p3);

	return BBox(min, max);
}

const Vector3<real> Parallelogram::ParallelogramValue::sample(real u1, real u2, Vector3<real> *N) const {
	Vector3<real> p(base + u1 * u + u2 * v);
	*N = uvw.w();

	return p;
}

const ONB &Parallelogram::ParallelogramValue::sample(real u1, real u2, Vector3<real> &pos) const {
	pos = base + u1 * u + u2 * v;
	return uvw;
}

const bool Parallelogram::ParallelogramValue::hit(const Ray& r, real &tHit, DifferentialGeometry& geometry) const {

	real tmax = r.tMax(), tmin = r.tMin();
	real dot1 = dot(r.direction(), norm);

	// check to see if ray is parallel to surface plane
	if(dot1 < Parallelogram::ParallelEpsilon() && dot1 > -Parallelogram::ParallelEpsilon()) return false;

	// find distance to surface plane
	real dot2 = dot(norm, base);
	real t = (dot2 - dot(norm, r.origin())) / dot1;
	if( t > tmax || t < tmin) return false;

	// see if point of plane intersection is within Parallelogram::ParallelogramValue
	Vector3<real> hit_plane(r.origin() + r.direction() * t);
	Vector3<real> offset(hit_plane - base);

	real u1 = dot(unorm, offset) / u.mag();
	if(u1 < -Parallelogram::Epsilon() || u1 > 1.0f+Parallelogram::Epsilon()) return false;

	real v1 = dot(vnorm, offset) / v.mag();
	if(v1 < -Parallelogram::Epsilon() || v1 > 1.0f+Parallelogram::Epsilon()) return false;

	// fill geometry record	
	geometry.hit_p = hit_plane;
	geometry.uvw = uvw;
	geometry.uv = v1 * uv2 + (1.0f - v1) * uv0 + u1 * uv1;
	geometry.shadingNormal = uvw.w();

	tHit = t;

	return true;
}

const bool Parallelogram::ParallelogramValue::shadowHit(const Ray& r) const {

	real tmin = r.tMin(), tmax = r.tMax();
	real dot1 = dot(r.direction(), norm);

	// check to see if ray is parallel to surface plane
	if(dot1 < Parallelogram::ParallelEpsilon() && dot1 > -Parallelogram::ParallelEpsilon()) return false;

	// find distance to surface plane
	real dot2 = dot(norm, base);
	real t = (dot2 - dot(norm, r.origin())) / dot1;
	if( t > tmax || t < tmin) return false;

	// see if point of plane intersection is within Parallelogram::ParallelogramValue
	Vector3<real> hit_plane(r.origin() + r.direction() * t);
	Vector3<real> offset(hit_plane - base);

	real u1 = dot(unorm, offset) / u.mag();
	if(u1 < -Parallelogram::Epsilon() || u1 > 1.0f+Parallelogram::Epsilon()) return false;

	real v1 = dot(vnorm, offset) / v.mag();
	return (v1 >= -Parallelogram::Epsilon() || v1 <= 1.0f+Parallelogram::Epsilon());
}
