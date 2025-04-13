#include <stdafx.h>
#include <Illumination/IlluminationHeader/Shape/Shape.h>

using namespace illuminationspace;

ShapeSet::ShapeSetValue::ShapeSetValue() : ShapeValue() {}
ShapeSet::ShapeSetValue::ShapeSetValue(const Transform &o2w) : ShapeValue(o2w) {}
ShapeSet::ShapeSetValue::ShapeSetValue(const std::vector<Shape> &s, const Transform &o2w)
: ShapeValue(o2w) {
	shapes = s;
	
	for (unsigned int i = 0; i < shapes.size(); ++i) 
		if(!shapes[i].canIntersect())
			basis::throwError("ShapeSet can include Intersectable Shapes only\n");

	calcAreas();
}

void ShapeSet::ShapeSetValue::calcAreas() {
	_area = 0;
	std::vector<real> areas;
	for (unsigned int i = 0; i < shapes.size(); ++i) {
		real a = shapes[i].area();
		_area += a;
		areas.push_back(a);
	}

	real prevCDF = 0;
	for (unsigned int i = 0; i < shapes.size(); ++i) {
		areaCDF.push_back(prevCDF + areas[i] / _area);
		prevCDF = areaCDF[i];
	}
}

const BBox ShapeSet::ShapeSetValue::objectBound() const {
	BBox ob;
	for (unsigned int i = 0; i < shapes.size(); ++i)
		ob.surround(shapes[i].objectBound());
	return ob;
}

const bool ShapeSet::ShapeSetValue::canIntersect() const {
	//for (unsigned int i = 0; i < shapes.size(); ++i)
	//	if (!shapes[i].canIntersect()) 
	//		return false;
	return true;
}

void ShapeSet::ShapeSetValue::refine(std::vector<Shape> &refined) const {
	for (unsigned int i = 0; i < shapes.size(); ++i) {
		if (shapes[i].canIntersect())	refined.push_back(shapes[i]);
		else							shapes[i].refine(refined);
	}
}

const Vector3<real> ShapeSet::ShapeSetValue::sample(real u1, real u2, Vector3<real> *Ns) const {
	real ls = random();
	unsigned int sn;
	for (sn = 0; sn < shapes.size()-1; ++sn)
		if (ls < areaCDF[sn]) 
			break;
	return shapes[sn].sample(u1, u2, Ns);
}

const bool ShapeSet::ShapeSetValue::hit(const Ray &ray, real &t_hit, DifferentialGeometry &dg) const {

	real old_tmax = ray.tMax();

	bool anyHit = false;
	real t;
	for (unsigned int i = 0; i < shapes.size(); ++i) {
		if (shapes[i].hit(ray, t, dg)) {
			anyHit = true;
			t_hit = t;
			ray.setTmax(t);
		}
	}

	ray.setTmax(old_tmax);

	return anyHit;
}

const bool ShapeSet::ShapeSetValue::shadowHit(const Ray &ray) const {
	bool anyHit = false;
	for (unsigned int i = 0; i < shapes.size(); ++i) {
		if (shapes[i].shadowHit(ray)) 
			return true;			
	}
	return false;
}	