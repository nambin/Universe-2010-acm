#include <stdafx.h>
#include <Principle/Geometry/Surface/OBJStruct.h>
#include <Illumination/IlluminationHeader/Shape/TriangleMesh.h>
#include <Illumination/IlluminationHeader/Shape/Triangle.h>

using namespace illuminationspace;

void TriangleMesh::loadFromOBJ(const char* const filename) {
	dynamic_cast<TriangleMeshValue*>(const_cast<ShapeValue*>(value.getptr()))->loadFromOBJ(filename);
}

TriangleMesh::TriangleMeshValue::TriangleMeshValue(const Transform &o2w, int nt, int nv, const int *vi, const Vector3<real> *P, const Vector3<real> *N, const Vector3<real> *S, const Vector2<real> *uv)
: Shape::ShapeValue(o2w) {

	numTriangles = nt;
	numVertices = nv;

	vertIdx.resize(numTriangles,3);
	vertIdx.set(vi);

	posVerts.resize(numVertices);
	for(int i=0;i<numVertices;++i)	
		posVerts[i] = ObjectToWorld.transformPoint(P[i]);

	uvw.resize(numTriangles);	
	for(int i=0;i<numTriangles;++i)	{
		Vector3<real> &p0 = posVerts[ vertIdx(i,0) ];
		Vector3<real> &p1 = posVerts[ vertIdx(i,1) ];
		Vector3<real> &p2 = posVerts[ vertIdx(i,2) ];
		uvw[i].initFromW(cross((p1-p0),(p2-p0)));
	}

	if(N != NULL)	{
		numNormals = numVertices;
		normals.resize(numNormals);
		normals.set(N);

		normalIdx.resize_set(vertIdx);
	}
	else 
		makeVertexNormals();

	if(S != NULL)	{
		s.resize(numVertices);
		s.set(S);
	}
	if(uv != NULL)	{
		uvs.resize(numVertices);
		uvs.set(uv);
	}
}

void TriangleMesh::TriangleMeshValue::loadFromOBJ(const char* const filename) {
	geometryspace::OBJStruct obj;
	obj.loadFromOBJ(filename);
	obj.unitize();

	if(obj.numVertices != obj.numTexCoords && obj.numTexCoords != 0)
		basis::throwError("NumVertices != NumNormals && obj.numNormals != 0 _ TriangleMesh::TriangleMeshValue::loadFromOBJ");

	numVertices = obj.numVertices;
	numTriangles = obj.numTriangles;
	numNormals = obj.numNormals;

	vertIdx.resize(numTriangles, 3);
	for(int i=0;i<numTriangles;++i) for(int j=0;j<3;++j) vertIdx(i,j) = obj.triangles[i].vIndex[j];
	
	posVerts.resize(obj.posVerts.size());
	for(int i=0;i<numVertices;++i) posVerts[i] = obj.posVerts[i];
	
	normalIdx.resize(numTriangles, 3);
	for(int i=0;i<numTriangles;++i) for(int j=0;j<3;++j) normalIdx(i,j) = obj.triangles[i].nIndex[j];

	normals.resize(obj.normals.size());
	for(int i=0;i<numNormals;++i) normals[i] = obj.normals[i];

	uvs.resize(obj.texCoords.size());
	for(unsigned int i=0;i<obj.texCoords.size();++i) uvs[i] = obj.texCoords[i];

	uvw.resize(numTriangles);	
	for(int i=0;i<numTriangles;++i)	{
		Vector3<real> &p0 = posVerts[ vertIdx(i,0) ];
		Vector3<real> &p1 = posVerts[ vertIdx(i,1) ];
		Vector3<real> &p2 = posVerts[ vertIdx(i,2) ];
		uvw[i].initFromW(cross((p1-p0),(p2-p0)));
	}
}

const BBox TriangleMesh::TriangleMeshValue::objectBound() const {
	BBox bbox;
	for(int i=0;i<numVertices;++i)
		bbox.surround(WorldToObject.transformPoint(posVerts[i]));
	return bbox;
}

const BBox TriangleMesh::TriangleMeshValue::worldBound() const {
	BBox bbox;
	for(int i=0;i<numVertices;++i)	
		bbox.surround(posVerts[i]);
	return bbox;
}

void TriangleMesh::TriangleMeshValue::refine(std::vector<Shape> &refined) const {
	for(int i=0;i<numTriangles;++i)	{
		Triangle t(ObjectToWorld, WorldToObject, this, i);
		refined.push_back(t);
	}
}

// Need to Modify
void TriangleMesh::TriangleMeshValue::makeVertexNormals() {
	numNormals = numVertices;
	normalIdx.resize_set(vertIdx);

	normals.resize(numVertices);

	std::vector<int> *t2n = new std::vector<int>[numVertices];
	for(int i=0;i<numTriangles;++i) {
		t2n[vertIdx(i,0)].push_back(i);
		t2n[vertIdx(i,1)].push_back(i);
		t2n[vertIdx(i,2)].push_back(i);
	}

	for(int i=0;i<numVertices;++i) {
		normals[i].set(0.0f);
		for(unsigned int j=0;j<t2n[i].size();++j) {
			normals[i] += uvw[t2n[i][j]].w();
		}
		normalize(normals[i]);
	}

	delete[] t2n;
}