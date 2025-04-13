#ifndef _ILLUMINATION_TRIANGLE_MESH_H_
#define _ILLUMINATION_TRIANGLE_MESH_H_

#include "../../../Basis/Basis/Array.h"
#include "../../../Basis/Basis/Array2.h"
#include "Shape.h"

namespace illuminationspace {

using basis::Array;
using basis::Array2;

class Triangle;

class ILLUMINATIONDLL TriangleMesh : public Shape {
public : 
	TriangleMesh() : Shape(new TriangleMeshValue()) {}
	TriangleMesh(const Transform &o2w, int nt, int nv, const int *vi, const Vector3<real> *P, 
				 const Vector3<real> *N=NULL, const Vector3<real> *S=NULL, const Vector2<real> *uv=NULL) 
	: Shape(new TriangleMeshValue(o2w, nt, nv, vi, P, N, S, uv)) {}

	void loadFromOBJ(const char* const filename);	

protected : 
	class ILLUMINATIONDLL TriangleMeshValue : public Shape::ShapeValue {
	public : 
		TriangleMeshValue() : Shape::ShapeValue() {}
		TriangleMeshValue(const Transform &o2w, int nt, int nv, const int *vi, const Vector3<real> *P, 
						  const Vector3<real> *N=NULL, const Vector3<real> *S=NULL, const Vector2<real> *uv=NULL);

		TriangleMeshValue* clone() const { return new TriangleMeshValue(*this); }

		void loadFromOBJ(const char* const filename);

		const bool canIntersect() const { return false; }

		const BBox objectBound() const;
		const BBox worldBound() const;
		
		void refine(std::vector<Shape> &refined) const;

		// Need to Modify
		void makeVertexNormals();

	public : 
		int numTriangles, numVertices, numNormals;

		Array2<int>				vertIdx;
		Array<Vector3<real>>	posVerts;

		Array<ONB>				uvw;

		Array2<int>				normalIdx;
		Array<Vector3<real>>	normals;

		Array<Vector2<real>>	uvs;
		Array<Vector3<real>>	s;
	};

public : 
	friend class Triangle;	
};

};		// illuminationspace

#endif