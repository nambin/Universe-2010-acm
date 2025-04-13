#include <stdafx.h>
#include <Principle/VirtualSpace/Scene/SceneLight.h>

using namespace virtualspace;

SceneLight_IOHandler::SceneLight_IOHandler(SceneLight *_s) {
	sceneLight = _s;
}

SceneLight_IOHandler::~SceneLight_IOHandler() {
}

void SceneLight_IOHandler::keyboardHandler(unsigned char key, int x, int y) {
}

void SceneLight_IOHandler::mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y) {
}

void SceneLight_IOHandler::mouseMotionHandler(int x, int y) {
}
