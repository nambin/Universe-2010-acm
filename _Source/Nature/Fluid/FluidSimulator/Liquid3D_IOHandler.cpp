#include <stdafx.h>
#include <Principle/VirtualSpace/Scene/Scene.h>
#include <Nature/Fluid/FluidSimulator/Liquid3D.h>

using namespace fluidspace;
using namespace virtualspace;

Liquid3D_IOHandler::Liquid3D_IOHandler(Liquid3D *_s) {
	sim = _s;
}

Liquid3D_IOHandler::~Liquid3D_IOHandler() {
}

void Liquid3D_IOHandler::keyboardHandler(unsigned char key, int x, int y) {
	switch(key) {
		case 'm' : sim->dMesh = !sim->dMesh; break;
		case 'b' : sim->dBoundingBox = !sim->dBoundingBox; break;
		case 's' : sim->dSurface = !sim->dSurface; break;
		case 'i' : sim->dInnerBoundary = !sim->dInnerBoundary; break;
		case 'x' : sim->tracking->exportOBJ(output_folder + "Liquid3D.obj"); cout << "Export OBJ Done" << endl; break;
		case ' ' : sim->isPause = !sim->isPause; break;

		case 'o' : scene->getCamera()->set(sim->viewEye(), sim->viewRef(), sim->viewUp(), sim->viewTheta()); break;		
	}
}

void Liquid3D_IOHandler::mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y) {
}

void Liquid3D_IOHandler::mouseMotionHandler(int x, int y) {
}
