#include <stdafx.h>
#include <Principle/Geometry/Surface/TriMeshConnectivity.h>

using namespace geometryspace;

TriMeshConnectivity::TriMeshConnectivity()	{ 
	init();				
}

TriMeshConnectivity::TriMeshConnectivity(const TriMeshConnectivity *mesh)	{ 
	init();	
	set(mesh);	
}

TriMeshConnectivity::TriMeshConnectivity(const TriMeshConnectivity &mesh)	{ 
	init();	
	set(mesh);	
}

void TriMeshConnectivity::init() {
	numVertices = numEdges = numFaces = 0;
	isBdryVertex.clear();
	v2f.clear(); e2f.clear(); f2f.clear(); 
	v2e.clear(); f2e.clear();;
	f2v.clear(); e2v.clear(); v2v.clear();
}

void TriMeshConnectivity::init(int nV, int nF, int nE) {
	numVertices = nV;
	numFaces = nF;
	numEdges = nE;	

	isBdryVertex.resize(numVertices);
	v2f.resize(numFaces,3);		e2f.resize(numFaces,3);		f2f.resize(numFaces,3);
	v2e.resize(numEdges,2);		f2e.resize(numEdges,2);
	f2v.resize(numVertices);	e2v.resize(numVertices);	v2v.resize(numVertices);	
}

void TriMeshConnectivity::set(const TriMeshConnectivity *mesh) {
	if(mesh == NULL) {
		release();
		return;
	}
	set(*mesh);
}

void TriMeshConnectivity::set(const TriMeshConnectivity &mesh) {
	release();	
	numVertices = mesh.numVertices;
	numFaces = mesh.numFaces;
	numEdges = mesh.numEdges;

	isBdryVertex.resize_set(mesh.isBdryVertex);
	v2f.resize_memcpy(mesh.v2f);	e2f.resize_memcpy(mesh.e2f);	f2f.resize_memcpy(mesh.f2f);
	v2e.resize_memcpy(mesh.v2e);	f2e.resize_memcpy(mesh.f2e);
	f2v.resize_set(mesh.f2v);		e2v.resize_set(mesh.e2v);		v2v.resize_set(mesh.v2v);
}

void TriMeshConnectivity::set(int nV, const std::vector<int> &_v2f) {
	release();

	numVertices = nV;
	numFaces = (int)((int)_v2f.size() / 3);

	v2f.resize(numFaces,3);
	v2f.set(_v2f);

	makeConnectivityData();
}

void TriMeshConnectivity::set(int nV, int nF, const int *_v2f) {
	release();

	numVertices = nV;
	numFaces = nF;

	v2f.resize(numFaces,3);
	v2f.set(_v2f);

	makeConnectivityData();
}

void TriMeshConnectivity::set(int nV, const Array2<int> &_v2f) {
	release();

	numVertices = nV;
	numFaces = static_cast<int>(_v2f.dimX());

	v2f.resize_memcpy(_v2f);

	makeConnectivityData();
}

void TriMeshConnectivity::release() {
	init();
}