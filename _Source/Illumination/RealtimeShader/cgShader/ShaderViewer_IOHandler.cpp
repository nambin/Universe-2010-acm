#include <stdafx.h>
#include <Illumination/RealtimeShader/ShaderViewer/ShaderViewer.h>

using namespace shaderspace;
using namespace virtualspace;

ShaderViewer_IOHandler::ShaderViewer_IOHandler(ShaderViewer *_s) {
	viewer = _s;
}

ShaderViewer_IOHandler::~ShaderViewer_IOHandler() {
}

void ShaderViewer_IOHandler::keyboardHandler(unsigned char key, int x, int y) {
	switch(key) {
		case 'w':
			viewer->vr.stepsize += 1.0/2048.0;
			if(viewer->vr.stepsize > 0.25) viewer->vr.stepsize = 0.25;
			break;
		case 'e':
			viewer->vr.stepsize -= 1.0/2048.0;
			if(viewer->vr.stepsize <= 1.0/200.0) viewer->vr.stepsize = 1.0/200.0;
			break;
		case 'o' : ;
			//scene->getCamera()->set(sim->viewEye(), sim->viewRef(), sim->viewUp(), sim->viewTheta()); break;
	}
}

void ShaderViewer_IOHandler::mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y) {
}

void ShaderViewer_IOHandler::mouseMotionHandler(int x, int y) {
}
