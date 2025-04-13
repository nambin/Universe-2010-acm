#include <stdafx.h>
#include <Principle/Virtualspace/Scene/Scene.h>
#include <Nature/Fluid/FluidScene/FluidScene.h>

using namespace fluidspace;
using namespace virtualspace;

FluidScene_IOHandler::FluidScene_IOHandler(FluidScene *_s) {
	sim = _s;
	selectedObject = NULL;
}

FluidScene_IOHandler::~FluidScene_IOHandler() {
}

void FluidScene_IOHandler::keyboardHandler(unsigned char key, int x, int y) {
	switch(key) {
		case 'u' : sim->updateOneStep(); break;
		case ' ' : sim->isPause = !sim->isPause; break;
		
		case 'o' : scene->getCamera()->set(sim->viewEye(), sim->viewRef(), sim->viewUp(), sim->viewTheta()); break;
		case '0' : scene->getCamera()->set(sim->viewEye_Side(), sim->viewRef_Side(), sim->viewUp_Side(), sim->viewTheta_Side()); break;
		
		case 'c' : scene->getCamera()->printCamera(); break;

		case '1' :	if(sim->fluid) selectedObject = dynamic_cast<SceneObject*>(sim->fluid);
					break;

		case '2' :	if(sim->objects) selectedObject = dynamic_cast<SceneObject*>(sim->objects);
					break;

		default  : if(selectedObject) selectedObject->getHandler()->keyboardHandler(key, x, y);
	}
}

void FluidScene_IOHandler::mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y) {
}

void FluidScene_IOHandler::mouseMotionHandler(int x, int y) {
}
