#include <stdafx.h>
#include <Principle/VirtualSpace/Camera/Camera.h>

using namespace virtualspace;

Camera_IOHandler::Camera_IOHandler(Camera *c) { camera = c; }
Camera_IOHandler::~Camera_IOHandler() {}

void Camera_IOHandler::keyboardHandler(unsigned char key, int x, int y) {
}

void Camera_IOHandler::mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y) {
	if(button == IO_MOUSE_RIGHT_BUTTON && state == IO_MOUSE_DOWN_STATE)		camera->trackballMouseDn(x, y);
	if(button == IO_MOUSE_RIGHT_BUTTON && state == IO_MOUSE_UP_STATE)		camera->trackballMouseUp(x, y);
	if(button == IO_MOUSE_MIDDLE_BUTTON && state == IO_MOUSE_DOWN_STATE)	camera->translateMouseDn(x, y);
	if(button == IO_MOUSE_MIDDLE_BUTTON && state == IO_MOUSE_UP_STATE)		camera->translateMouseUp(x, y);
	if(button == IO_MOUSE_LEFT_BUTTON && state == IO_MOUSE_DOWN_STATE)		camera->panMouseDn(x, y);
	if(button == IO_MOUSE_LEFT_BUTTON && state == IO_MOUSE_UP_STATE)		camera->panMouseUp(x, y);
	if(button == IO_MOUSE_WHEEL_UP)											camera->panin();
	if(button == IO_MOUSE_WHEEL_DOWN)										camera->panout();
}

void Camera_IOHandler::mouseMotionHandler(int x, int y) {
	camera->mouseMove(x, y);
}