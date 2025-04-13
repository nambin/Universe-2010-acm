#ifndef _ILLUMINATION_TRIANGLE_H_
#define _ILLUMINATION_TRIANGLE_H_

#include "Shape.h"
#include "TriangleMesh.h"

namespace illuminationspace {

class ILLUMINATIONDLL Triangle : public Shape {
public :
	Triangle(const Transform &o2w, const Transform &w2o, const TriangleMesh::TriangleMeshValue *m, int i)
	: Shape(new TriangleValue(o2w, w2o, m, i)) {}

protected : 
	class ILLUMINATIONDLL TriangleValue : public Shape::ShapeValue {
	public :
		TriangleValue(const Transform &o2w, const Transform &w2o, const TriangleMesh::TriangleMeshValue *m, int i)
		: meshptr(const_cast<TriangleMesh::TriangleMeshValue*>(m)), v(&meshptr->vertIdx.getptr()[i*3]), n(&meshptr->normalIdx.getptr()[i*3]), uvw(&meshptr->uvw.getptr()[i]) {
			ObjectToWorld = o2w;
			WorldToObject = w2o; 
		}

		TriangleValue* clone() const { return new TriangleValue(*this); }

		const BBox objectBound() const;
		const BBox worldBound() const;

		const real area() const;

		const Vector3<real> sample(real u1, real u2, Vector3<real> *N) const;

		const bool hit(const Ray& r, real &tHit, DifferentialGeometry &geometry) const;
		const bool shadowHit(const Ray& r) const;

	public :		
		const basis::RCPtr<TriangleMesh::TriangleMeshValue> meshptr;
		const int *v;
		const int *n;
		const ONB *uvw;
	};
};

};		// illuminationspace

#endif