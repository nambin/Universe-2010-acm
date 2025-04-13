#include <stdafx.h>
#include <Illumination/IlluminationHeader/Shape/Shape.h>

#ifdef _WIN32
#pragma warning(disable:4723)   
#endif

using namespace std;
using namespace illuminationspace;

Shape::ShapeValue::ShapeValue(const Transform &o2w) {
	ObjectToWorld.set(o2w.getMatrixRef(), o2w.getMatrixInvRef());
	WorldToObject.set(o2w.getMatrixInvRef(), o2w.getMatrixRef());
}	

const BBox Shape::ShapeValue::worldBound() const { 
	return ObjectToWorld.transformBBox(objectBound()); 
}

void Shape::ShapeValue::refine(std::vector<Shape> &refined) const {	
	basis::throwError("Not Implement Shape::ShapeValue::refine method\n"); 
}

void Shape::ShapeValue::getShadingGeometry(const Transform &obj2world, const DifferentialGeometry &dg, DifferentialGeometry &dgShading) const {
	dgShading = dg;
}

const real Shape::ShapeValue::area() const { 
	basis::throwError("Not Implement Shape::ShapeValue::area method\n");
	return 0.0f; 
}

const real Shape::ShapeValue::pdf() const {
	return 1.0f / area();
} 

const Vector3<real> Shape::ShapeValue::sample(real u1, real u2, Vector3<real> *N) const { 
	basis::throwError("Not Implement Shape::ShapeValue::sample method\n");
	return Vector3<real>(0,0,0); 
}

const Vector3<real> Shape::ShapeValue::sample(const Vector3<real> &pt, real u1, real u2, Vector3<real> *N) const {
	return sample(u1, u2, N);
}

const ONB &Shape::ShapeValue::sample(real u1, real u2, Vector3<real> &pos)const {
	basis::throwError("Not Implement Shape::ShapeValue::sample(u1,u2,pos,uvw) method\n");
	exit(0);
}

const bool Shape::ShapeValue::hit(const Ray& r, real& tHit, DifferentialGeometry& geometry) const { 
	basis::throwError("Not Implement Shape::ShapeValue::hit method\n");
	return false; 
}

const bool Shape::ShapeValue::shadowHit(const Ray& r) const { 
	basis::throwError("Not Implement Shape::ShapeValue::shadowHit method\n"); 
	return false; 
}