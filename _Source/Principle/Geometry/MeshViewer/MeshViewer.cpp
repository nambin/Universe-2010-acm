#include <stdafx.h>
#include <Principle/Geometry/MeshViewer/MeshViewer.h>

using namespace geometryspace;

MeshViewer::MeshViewer() {
	ioHandler = new MeshViewer_IOHandler(this);
	selected_face= 800;
}

MeshViewer::MeshViewer(const TriMesh<Vector3f> &_mesh) {
	set(_mesh);
	ioHandler = new MeshViewer_IOHandler(this);
	selected_face = 800;
}

MeshViewer::~MeshViewer() {
}

void MeshViewer::set(const TriMesh<Vector3f> &_mesh) {
	mesh.set(_mesh);
}

void MeshViewer::load(const char * const filename, bool change_order) {
	
	mesh_simple.load(filename);
	mesh_simple.unitize();
	if(change_order)
		mesh_simple.changeTriangleOrder();
	mesh_simple.makeFaceNormals();
	mesh_simple.makeVertexNormals();

	//mesh.load(filename);
	//mesh.unitize();	
	//if(change_order)
	//	mesh.changeTriangleOrder();
	//mesh.makeVertexNormals();	
		
	//mesh.makeConnectivityData();
	//std::cout << mesh.isClosedTwoManifold() << std::endl;
}

void MeshViewer::glBind() const {

	//SceneObject::glBind();
	mesh_simple.glDraw_SmoothSurface(RGBA::white() * 0.7f);

	//mesh.glDrawTriangle(selected_face);
	//mesh.glDrawVertex(mesh.v2f(selected_face, 0), 5.0f, RGBA::red());
	//mesh.glDrawVertex(mesh.v2f(selected_face, 1), 5.0f, RGBA::green());
	//mesh.glDrawVertex(mesh.v2f(selected_face, 2), 5.0f, RGBA::blue());
	
	mesh.glDrawSmoothSurface(RGBA::white() * 0.7f);	
	//mesh.glDrawLineStrip();
	//mesh.glDrawNormalMap();
	//mesh.glDraw();	

	//objStruct.glDraw();

	//plyStruct.glDrawSimpleMesh();
	//plyStruct.glDrawVertexColoredMesh();
}

void MeshViewer::idle() {
}