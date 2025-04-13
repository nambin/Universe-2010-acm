#include <stdafx.h>
#include <Principle/VirtualSpace/Scene/Scene.h>
#include <Nature/Fluid/FluidSimulator/Liquid2D.h>

using namespace fluidspace;
using namespace virtualspace;

Liquid2D_IOHandler::Liquid2D_IOHandler(Liquid2D *_s) {
	sim = _s;
}

Liquid2D_IOHandler::~Liquid2D_IOHandler() {
}

void Liquid2D_IOHandler::keyboardHandler(unsigned char key, int x, int y) {
	switch(key) {
		case 'b' : sim->dBoundingBox = !sim->dBoundingBox; break;
		case 'i' : sim->dInternalBoundary = !sim->dInternalBoundary; break;
		case 'g' : sim->dGrid = !sim->dGrid; break;
		case 'l' : sim->dLiquid = !sim->dLiquid; break;
		case 'f' : sim->dFilledInterface = !sim->dFilledInterface; break;
		case 'c' : sim->dParticle = !sim->dParticle; break;
		case 'p' : sim->dPhi = !sim->dPhi; break;
		case 'r' : sim->dPressure = !sim->dPressure; break;
		case 'v' : sim->dVelocity = !sim->dVelocity; break;
		case 'a' : sim->dAuxiliary = !sim->dAuxiliary; break;
		case 'u' : sim->updateOneStep(); break;
		case ' ' : sim->isPause = !sim->isPause; break;
		case '+' : sim->scaleVelocity *= 1.1; break;
		case '-' : sim->scaleVelocity /= 1.1; break;

		case 'o' : scene->getCamera()->set(sim->viewEye(), sim->viewRef(), sim->viewUp(), sim->viewTheta()); break;
	}
}

void Liquid2D_IOHandler::mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y) {
}

void Liquid2D_IOHandler::mouseMotionHandler(int x, int y) {
}
