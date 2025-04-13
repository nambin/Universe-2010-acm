#ifndef _VIRTUALSPACE_IOHANDLER_H_
#define _VIRTUALSPACE_IOHANDLER_H_

#include <iostream>
#include "../../Principle.h"

namespace virtualspace {

class Scene;
class PRINCIPLEDLL IOHandler {
public : 
	enum IO_MOUSE_BUTTON { IO_MOUSE_LEFT_BUTTON, IO_MOUSE_RIGHT_BUTTON, IO_MOUSE_MIDDLE_BUTTON, IO_MOUSE_WHEEL_UP, IO_MOUSE_WHEEL_DOWN };
	enum IO_MOUSE_STATE	 { IO_MOUSE_DOWN_STATE, IO_MOUSE_UP_STATE };

public : 
	IOHandler() {};
	virtual ~IOHandler();

	virtual void keyboardHandler(unsigned char key, int	x, int y) = 0;
	virtual void mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y) = 0;
	virtual void mouseMotionHandler(int x, int y) = 0;

	void setScene(const Scene *s) { scene = s; }
	const Scene* const getScene() const { return scene; }

protected : 
	const Scene *scene;
};

};		// namespace virtualspace

#endif