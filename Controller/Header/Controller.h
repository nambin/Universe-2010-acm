#ifndef _CONTROLLER_H_
#define _CONTROLLER_H_

#include <Principle/VirtualSpace/Scene/Scene.h>
#include <Basis/Timer/Timer.h>

namespace controlspace {

class Controller {
public :
	//////////////////////////////////////////////////////////////////////////
	// Constructor
	Controller();
	~Controller();

	void init();
	void set();

	void release();
	
	//////////////////////////////////////////////////////////////////////////
	// Render
	inline void render() const;	

	//////////////////////////////////////////////////////////////////////////
	// IdleFunction
	inline void idleFunction();

	//////////////////////////////////////////////////////////////////////////
	// Handler
	inline void keyboardHandler(unsigned char key, int	x, int y);
	void mouseHandler(int button, int state, int x, int y);
	inline void mouseMotionHandler(int x, int y);

	//////////////////////////////////////////////////////////////////////////
	// View Change
	inline void changeSize(unsigned int _width, unsigned int _height);	

private : 
	virtualspace::Scene *scene;
};

inline void Controller::render() const {	
	scene->glDraw();		
	glutSwapBuffers();
}

inline void Controller::idleFunction() {
	scene->idle();	
}

inline void Controller::keyboardHandler(unsigned char key, int x, int y) {
	scene->keyboardHandler(key, x, y);
}

inline void Controller::mouseMotionHandler(int x, int y) {
	scene->mouseMotionHandler(x,y);
}

inline void Controller::changeSize(unsigned int _width, unsigned int _height) {
	scene->resizeScreen(_width, _height);
}

};		// namespace controlspace;

#endif