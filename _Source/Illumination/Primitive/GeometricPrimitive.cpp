#include <stdafx.h>
#include <Illumination/IlluminationHeader/Primitive/GeometricPrimitive.h>

using namespace illuminationspace;

void GeometricPrimitive::GeometricPrimitiveValue::refine(std::vector<Primitive> &refined) const {
	std::vector<Shape> vecShape;
	shape.refine(vecShape);

	for(unsigned int i=0;i<vecShape.size();++i) {
		GeometricPrimitive prim(vecShape[i], material, areaLight);
		refined.push_back(prim);
	}
}

const bool GeometricPrimitive::GeometricPrimitiveValue::hit(const Ray& r, Intersection &intersection) const {
	if(!shape.hit(r, intersection.tHit, intersection.dg))
		return false;

	intersection.primitiveValue = this;
	return true;
}

const bool GeometricPrimitive::GeometricPrimitiveValue::shadowHit(const Ray& r) const {
	return shape.shadowHit(r);
}

