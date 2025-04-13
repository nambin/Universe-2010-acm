#include <stdafx.h>
#include <Principle/VirtualSpace/Scene/Scene.h>
#include <Principle/Geometry/BiTree/OctreeViewer.h>

using namespace virtualspace;
using namespace geometryspace;

OctreeViewer_IOHandler::OctreeViewer_IOHandler(OctreeViewer *_viewer) {
	viewer = _viewer;
}

OctreeViewer_IOHandler::~OctreeViewer_IOHandler() {}

void OctreeViewer_IOHandler::keyboardHandler(unsigned char key, int x, int y) {	

	switch(key) {
		case 'o' : scene->getCamera()->set(viewer->viewEye(), viewer->viewRef(), viewer->viewUp(), viewer->viewTheta()); break;
		case '+' : if(viewer->idxCell < viewer->octree.cells.size()-1) viewer->idxCell++; break;
		case '-' : if(viewer->idxCell > 0) viewer->idxCell--; break;
		case 't' : viewer->drawTree = !viewer->drawTree; break;
		case 'l' : viewer->drawOnlyMinimumLeaf = !viewer->drawOnlyMinimumLeaf; break;
		case 'r' : viewer->drawRenderingPrimitive = !viewer->drawRenderingPrimitive; break;
	}
}

void OctreeViewer_IOHandler::mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y) {
}

void OctreeViewer_IOHandler::mouseMotionHandler(int x, int y) {
}