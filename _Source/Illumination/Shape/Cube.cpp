#include <stdafx.h>
#include <Illumination/IlluminationHeader/Shape/Cube.h>
#include <Illumination/IlluminationHeader/Shape/Pallelogram.h>

using namespace illuminationspace;

Cube::CubeValue::CubeValue(const Transform &o2w, const Vector3<real>& v0, const Vector3<real>& v1)
: ShapeSet::ShapeSetValue(o2w) {
	Transform identity = o2w;

	Vector3<real> _min = minVec(v0, v1);
	Vector3<real> _max = maxVec(v0, v1);	

	Vector3<real> diff = _max - _min;
	Vector3<real> x(diff[0], 0, 0);
	Vector3<real> y(0, diff[1], 0);
	Vector3<real> z(0, 0, diff[2]);

	shapes.push_back(Parallelogram(identity, _min  , y, x));	
	shapes.push_back(Parallelogram(identity, _min+z, x, y));	
	shapes.push_back(Parallelogram(identity, _min  , z, y));	
	shapes.push_back(Parallelogram(identity, _min+x, y, z));	
	shapes.push_back(Parallelogram(identity, _min  , x, z));	
	shapes.push_back(Parallelogram(identity, _min+y, z, x));	

	calcAreas();

	cube.set(_min, _max);
}

const BBox Cube::CubeValue::objectBound() const {	
	return cube;
}

const bool Cube::CubeValue::canIntersect() const {	
	return true;
}

void Cube::CubeValue::refine(std::vector<Shape> &refined) const {
	Shape::ShapeValue::refine(refined);
}

const Vector3<real> Cube::CubeValue::sample(real u1, real u2, Vector3<real> *Ns) const {
	real ls = random();
	unsigned int sn;
	for (sn = 0; sn < shapes.size()-1; ++sn)
		if (ls < areaCDF[sn]) 
			break;
	return shapes[sn].sample(u1, u2, Ns);
}

const bool Cube::CubeValue::hit(const Ray &ray, real &t_hit, DifferentialGeometry &dg) const {

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

const bool Cube::CubeValue::shadowHit(const Ray &ray) const {
	bool anyHit = false;
	for (unsigned int i = 0; i < shapes.size(); ++i) {
		if (shapes[i].shadowHit(ray)) 
			return true;			
	}
	return false;
}	