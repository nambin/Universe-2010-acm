#include <stdafx.h>
#include <Principle/VirtualSpace/Scene/Scene.h>
#include <Principle/Geometry/BiTree/QuadTreeViewer.h>

using namespace virtualspace;
using namespace geometryspace;

QuadTreeViewer_IOHandler::QuadTreeViewer_IOHandler(QuadTreeViewer *_viewer) {
	viewer = _viewer;
}

QuadTreeViewer_IOHandler::~QuadTreeViewer_IOHandler() {}

void QuadTreeViewer_IOHandler::keyboardHandler(unsigned char key, int x, int y) {	

	switch(key) {
		case 'o' : scene->getCamera()->set(viewer->viewEye(), viewer->viewRef(), viewer->viewUp(), viewer->viewTheta()); break;
		case '+' : if(viewer->idxCell < viewer->quadtree.cells.size()-1) viewer->idxCell++; break;
		case '-' : if(viewer->idxCell > 0) viewer->idxCell--; break;
		case 't' : viewer->drawTree = !viewer->drawTree; break;
		case 'r' : viewer->drawRenderingPrimitive = !viewer->drawRenderingPrimitive; break;
	}
}

void QuadTreeViewer_IOHandler::mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y) {
}

void QuadTreeViewer_IOHandler::mouseMotionHandler(int x, int y) {
}