#include <stdafx.h>
#include <Principle/VirtualSpace/Scene/Scene.h>
#include <Principle/Geometry/MeshViewer/TriangleMeshViewer.h>

using namespace virtualspace;
using namespace geometryspace;

TriangleMeshViewer_IOHandler::TriangleMeshViewer_IOHandler(TriangleMeshViewer *_viewer) {
	viewer = _viewer;
}

TriangleMeshViewer_IOHandler::~TriangleMeshViewer_IOHandler() {}

void TriangleMeshViewer_IOHandler::keyboardHandler(unsigned char key, int x, int y) {	

	switch(key) {
		case 'o' : scene->getCamera()->set(viewer->viewEye(), viewer->viewRef(), viewer->viewUp(), viewer->viewTheta()); break;
		case 'i' : viewer->iv++; if(viewer->iv >= (unsigned int)viewer->mesh.numVertices) viewer->iv = 0; break;
		case 'u' : viewer->iv--; if(viewer->iv == -1) viewer->iv = viewer->mesh.numVertices-1; break;
		case 'n' : viewer->nth++;
	}
}

void TriangleMeshViewer_IOHandler::mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y) {
}

void TriangleMeshViewer_IOHandler::mouseMotionHandler(int x, int y) {
}