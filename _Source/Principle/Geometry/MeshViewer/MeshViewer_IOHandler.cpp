#include <stdafx.h>
#include <Principle/VirtualSpace/Scene/Scene.h>
#include <Principle/Geometry/MeshViewer/MeshViewer.h>

using namespace virtualspace;
using namespace geometryspace;

MeshViewer_IOHandler::MeshViewer_IOHandler(MeshViewer *_viewer) {
	viewer = _viewer;
}

MeshViewer_IOHandler::~MeshViewer_IOHandler() {}

void MeshViewer_IOHandler::keyboardHandler(unsigned char key, int x, int y) {	

	switch(key) {
		case 'o' : scene->getCamera()->set(viewer->viewEye(), viewer->viewRef(), viewer->viewUp(), viewer->viewTheta()); break;
		case '+' : viewer->selected_face++; cout << viewer->selected_face << endl; break;
		case '-' : viewer->selected_face--; cout << viewer->selected_face << endl; break;
	}
}

void MeshViewer_IOHandler::mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y) {
}

void MeshViewer_IOHandler::mouseMotionHandler(int x, int y) {
}