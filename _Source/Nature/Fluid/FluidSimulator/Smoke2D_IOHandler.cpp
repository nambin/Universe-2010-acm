#include <stdafx.h>
#include <Principle/VirtualSpace/Scene/Scene.h>
#include <Nature/Fluid/FluidSimulator/Smoke2D.h>

using namespace fluidspace;
using namespace virtualspace;

Smoke2D_IOHandler::Smoke2D_IOHandler(Smoke2D *_s) {
	sim = _s;
}

Smoke2D_IOHandler::~Smoke2D_IOHandler() {
}

void Smoke2D_IOHandler::keyboardHandler(unsigned char key, int x, int y) {
	switch(key) {
		case 'k' : sim->dKinematicBoundary = !sim->dKinematicBoundary; break;
		case 'b' : sim->dBoundingBox = !sim->dBoundingBox; break;
		case 'd' : sim->dDensity = !sim->dDensity; break;
		case 'p' : sim->dPressure = !sim->dPressure; break;
		case 'v' : sim->dVelocity = !sim->dVelocity; break;
		case ' ' : sim->isPause = !sim->isPause; break;
		case '+' : sim->scaleVelocity *= 1.1f; break;
		case '-' : sim->scaleVelocity /= 1.1f; break;

		case 'o' : scene->getCamera()->set(sim->viewEye(), sim->viewRef(), sim->viewUp(), sim->viewTheta()); break;
	}
}

void Smoke2D_IOHandler::mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y) {
}

void Smoke2D_IOHandler::mouseMotionHandler(int x, int y) {
}
