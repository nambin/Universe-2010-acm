#include <stdafx.h>
#include "Header/Controller.h"

using namespace virtualspace;
using namespace controlspace;

void Controller::mouseHandler(int button, int state, int x, int y) {

	IOHandler::IO_MOUSE_BUTTON _button;
	IOHandler::IO_MOUSE_STATE _state;

	switch(button) {
		case GLUT_RIGHT_BUTTON	: _button = IOHandler::IO_MOUSE_RIGHT_BUTTON;	break;
		case GLUT_LEFT_BUTTON	: _button = IOHandler::IO_MOUSE_LEFT_BUTTON;	break;
		case GLUT_MIDDLE_BUTTON	: _button = IOHandler::IO_MOUSE_MIDDLE_BUTTON;	break;
//		case GLUT_WHEEL_UP		: _button = IOHandler::IO_MOUSE_WHEEL_UP;		break;
//		case GLUT_WHEEL_DOWN	: _button = IOHandler::IO_MOUSE_WHEEL_DOWN;		break;
		default : 
			std::cerr << "Unknown Mouse Button !!" << std::endl;
			return;		
	}
	switch(state) {
		case GLUT_DOWN	: _state = IOHandler::IO_MOUSE_DOWN_STATE;	break;
		case GLUT_UP	: _state = IOHandler::IO_MOUSE_UP_STATE;	break;
		default : 
			std::cerr << "Unknown Mouse State !!" << std::endl;
			return;
	}

	scene->mouseHandler(_button, _state, x, y);
}