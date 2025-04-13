#include <stdafx.h>
#include <Principle/Geometry/MeshViewer/TriangleMeshViewer.h>

using namespace geometryspace;

TriangleMeshViewer::TriangleMeshViewer() {
	ioHandler = new TriangleMeshViewer_IOHandler(this);
}

TriangleMeshViewer::TriangleMeshViewer(const TriMesh<Vector2f> &_mesh) {
	set(_mesh);
	ioHandler = new TriangleMeshViewer_IOHandler(this);
}

TriangleMeshViewer::~TriangleMeshViewer() {
}

void TriangleMeshViewer::set(const TriMesh<Vector2f> &_mesh) {
	mesh.set(_mesh);
	mesh.makeConnectivityData();
	iv = nth = 0;
}

void TriangleMeshViewer::glBind() const {
	SceneObject::glBind();
	mesh.glDrawLineStrip2D(2);

	if(nth >= mesh.v2v[iv].size() && nth >= mesh.e2v[iv].size() && nth >= mesh.f2v[iv].size()) {
		cout << mesh.v2v[iv].size() << endl;
		cout << mesh.e2v[iv].size() << endl;
		cout << mesh.f2v[iv].size() << endl;
		nth = 0;
	}

	mesh.glDrawV2V_2D(iv, nth, 10.0f, RGBA::red(), RGBA::magenta());
	mesh.glDrawE2V_2D(iv, nth, 10.0f, RGBA::red(), RGBA::yellow());
	mesh.glDrawF2V_2D(iv, nth, 10.0f, RGBA::red(), RGBA::white());
//	mesh.glDrawBoundaryVertex2D(10.0f, RGBA::red());
}

void TriangleMeshViewer::idle() {
}
