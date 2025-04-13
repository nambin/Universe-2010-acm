#include <stdafx.h>
#include <Principle/Virtualspace/Scene/Scene.h>
#include <Nature/MassSpring/ClothScene/ClothScene.h>

using namespace springspace;
using namespace virtualspace;

ClothScene_IOHandler::ClothScene_IOHandler(ClothScene *_s) {
	sim = _s;
}

ClothScene_IOHandler::~ClothScene_IOHandler() {
}

void ClothScene_IOHandler::keyboardHandler(unsigned char key, int x, int y) {
	switch(key) {
		case 'u' : sim->updateOneStep(); break;
		case ' ' : sim->isPause = !sim->isPause; break;

		case 'f' : sim->dSurface = !sim->dSurface; break;
		case 't' : sim->dStretchSpring = !sim->dStretchSpring; break;
		case 's' : sim->dShearSpring = !sim->dShearSpring; break;
		case 'b' : sim->dBendingSpring = !sim->dBendingSpring; break;
		case 'p' : sim->dParticle = !sim->dParticle; break;
		case 'j' : sim->dObjects = !sim->dObjects; break;

		case 'o' : scene->getCamera()->set(sim->viewEye(), sim->viewRef(), sim->viewUp(), sim->viewTheta()); break;
	}
}

void ClothScene_IOHandler::mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y) {
}

void ClothScene_IOHandler::mouseMotionHandler(int x, int y) {
}
