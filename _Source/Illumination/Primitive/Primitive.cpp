#include <stdafx.h>
#include <Illumination/IlluminationHeader/Primitive/Primitive.h>

using namespace illuminationspace;

Primitive::PrimitiveValue* Primitive::PrimitiveValue::clone() const { 
	basis::throwError("Not Implement Primitive::PrimitiveValue::clone method\n"); 
	return NULL; 
}

void Primitive::PrimitiveValue::refine(std::vector<Primitive> &refined) const {
	throwError("Not Implement Primitive::PrimitiveValue::refine method\n"); 
}

void Primitive::PrimitiveValue::fullyRefine(std::vector<Primitive> &refined) const {

	Primitive prim0(const_cast<Primitive::PrimitiveValue*>(this));

	std::vector<Primitive> todo;
	todo.push_back(prim0);

	while(todo.size()) {
		Primitive prim = todo.back();
		todo.pop_back();

		if(prim.canIntersect()) refined.push_back(prim);
		else					prim.refine(todo);
	}
}

#pragma warning(disable:4172)   

const AreaLight& Primitive::PrimitiveValue::getAreaLight() const {
	basis::throwError("Not Implement Primitive::PrimitiveValue::getAreaLight method\n");  
	return AreaLight(); 
}

BSDF *Primitive::PrimitiveValue::getBSDF(const DifferentialGeometry &dg) const {
	basis::throwError("Not Implement Primitive::PrimitiveValue::getBSDF method\n");  
	return NULL; 
}

const Spectrum Primitive::PrimitiveValue::extinction(const DifferentialGeometry& dg, const Vector3<real>& wi) const {
	basis::throwError("Not Implement Primitive::PrimitiveValue::extinction method\n");  
	return Spectrum::black();
}

const bool Primitive::PrimitiveValue::hit(const Ray& r, Intersection &intersection) const {
	basis::throwError("Not Implement Primitive::PrimitiveValue::hit method\n");  
	return false; 
}

const bool Primitive::PrimitiveValue::shadowHit(const Ray& r) const {
	basis::throwError("Not Implement Primitive::PrimitiveValue::shadowHit method\n");
	return false; 
}
