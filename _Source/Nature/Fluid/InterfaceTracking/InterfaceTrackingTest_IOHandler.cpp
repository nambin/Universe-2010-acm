#include <stdafx.h>
#include <Principle/VirtualSpace/Scene/Scene.h>
#include <Nature/Fluid/InterfaceTracking/InterfaceTrackingTest.h>

using namespace fluidspace;
using namespace virtualspace;

InterfaceTrackingTest_IOHandler::InterfaceTrackingTest_IOHandler(InterfaceTrackingTest *_s) {
	sim = _s;
}

InterfaceTrackingTest_IOHandler::~InterfaceTrackingTest_IOHandler() {
}

void InterfaceTrackingTest_IOHandler::keyboardHandler(unsigned char key, int x, int y) {
	switch(key) {
		case 'b' : sim->dBoundingBox = !sim->dBoundingBox; break;
		case 'g' : sim->dGrid = !sim->dGrid; break;
		case 'l' : sim->dLiquid = !sim->dLiquid; break;		
		case 'p' : sim->dPhi = !sim->dPhi; break;
		case 'a' : sim->dAuxiliary = !sim->dAuxiliary; break;
		case 'f' : sim->dFilledMesh = !sim->dFilledMesh; break;
		case 'u' : sim->updateOneStep(); break;		
		case ' ' : sim->isPause = !sim->isPause; break;

		case 'o' : scene->getCamera()->set(sim->viewEye(), sim->viewRef(), sim->viewUp(), sim->viewTheta()); break;
	}
}

void InterfaceTrackingTest_IOHandler::mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y) {
}

void InterfaceTrackingTest_IOHandler::mouseMotionHandler(int x, int y) {
}
