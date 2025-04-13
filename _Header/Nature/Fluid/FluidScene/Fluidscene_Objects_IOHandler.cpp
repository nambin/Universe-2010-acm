#include <stdafx.h>
#include <Principle/Virtualspace/Scene/Scene.h>
#include <Nature/Fluid/FluidScene/FluidScene_Objects.h>

using namespace fluidspace;
using namespace virtualspace;

FluidScene_Objects_IOHandler::FluidScene_Objects_IOHandler(FluidScene_Objects *_s) {
	objs = _s;
}

FluidScene_Objects_IOHandler::~FluidScene_Objects_IOHandler() {
}

void FluidScene_Objects_IOHandler::keyboardHandler(unsigned char key, int x, int y) {
	switch(key) {
		case 'd' :	if(objs) objs->dObject = !objs->dObject;
					break;
	}
}

void FluidScene_Objects_IOHandler::mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y) {
}

void FluidScene_Objects_IOHandler::mouseMotionHandler(int x, int y) {
}
