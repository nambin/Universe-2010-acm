#include <stdafx.h>
#include <Illumination/IlluminationHeader/Shape/Sphere.h>

using namespace illuminationspace;

const BBox Sphere::SphereValue::objectBound() const {
	return BBox(center - Vector3<real>(radius), center + Vector3<real>(radius));
}

const bool Sphere::SphereValue::hit(const Ray& ray, real &tHit, DifferentialGeometry& geometry) const {

	Ray r = WorldToObject.transformRay(ray);
	//Ray r = ray;
	real tmin = ray.tMin(), tmax = ray.tMax();

	Vector3<real> tmp = r.origin() - center;

	double a = dot(r.direction(), r.direction());
	double b = 2.0f * dot(r.direction(), tmp);
	double c = dot(tmp, tmp) - radius * radius;

	real discriminant = static_cast<real>(b*b - 4.0f*a*c);

	// first check to see if ray intersects sphere
	if(discriminant > 0) {
		discriminant = sqrt(discriminant);
		double t = (-b - discriminant) / (2.0f*a);

		// now check for valid interval
		if(t < tmin)				t = (-b + discriminant) / (2.0f*a);
		if(t < tmin || t > tmax)	return false;

		// we have a valid hit
		tHit = static_cast<real>(t);
		geometry.hit_p = ray.origin() + tHit * ray.direction();

		Vector3<real> oNormal = unitVector(r.origin() + tHit * r.direction() - center);
		geometry.shadingNormal = ObjectToWorld.transformNormal(oNormal);
		geometry.uvw.initFromUnitVector(geometry.shadingNormal);
		geometry.uv = Vector2<real>(0,0);

		return true;
	}

	return false;
}

const bool Sphere::SphereValue::shadowHit(const Ray& ray) const {

	Ray r = WorldToObject.transformRay(ray);
	//Ray r = ray;
	real tmin = r.tMin(), tmax = r.tMax();

	Vector3<real> tmp = r.origin() - center;

	double a = dot(r.direction(), r.direction());
	double b = 2.0f * dot(r.direction(), tmp);
	double c = dot(tmp, tmp) - radius * radius;

	real discriminant = static_cast<real>(b*b - 4.0f*a*c);

	// first check to see if ray intersects sphere
	if(discriminant > 0) {
		discriminant = sqrt(discriminant);
		double t = (-b - discriminant) / (2.0f*a);

		// now check for valid interval
		if(t < tmin)				t = (-b + discriminant) / (2.0f*a);
		if(t < tmin || t > tmax)	return false;

		// we have a valid hit
		return true;
	}
	return false;
}
