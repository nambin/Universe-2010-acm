#include <stdafx.h>
#include <Principle/VirtualSpace/IOHandler/IOHandler.h>

using namespace virtualspace;

IOHandler::~IOHandler() {
}

void IOHandler::keyboardHandler(unsigned char key, int	x, int y) {
}

void IOHandler::mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y) {
}

void IOHandler::mouseMotionHandler(int x, int y) {
}