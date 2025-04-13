#ifndef _TEST_OPENGL_SCENE_H_
#define _TEST_OPENGL_SCENE_H_

#include <Principle/VirtualSpace/Scene/SceneObject.h>
#include <Principle/VirtualSpace/IOHandler/IOHandler.h>

using namespace basis;
using namespace virtualspace;

class TestOpenGLScene;
class TestOpenGLScene_IOHandler : public IOHandler {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	TestOpenGLScene_IOHandler(TestOpenGLScene *_sim);
	~TestOpenGLScene_IOHandler();

	//////////////////////////////////////////////////////////////////////////
	// Handler
	void keyboardHandler(unsigned char key, int	x, int y);
	void mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y);
	void mouseMotionHandler(int x, int y);

private : 
	TestOpenGLScene*		sim;	
};

class TestOpenGLScene : public SceneObject {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	TestOpenGLScene() { ioHandler = new TestOpenGLScene_IOHandler(this); }
	~TestOpenGLScene() {}

	//////////////////////////////////////////////////////////////////////////
	// glBind
	virtual void glBind() const;
	void glDraw() const;

	void resizeScreen(unsigned int width, unsigned int height);

	//////////////////////////////////////////////////////////////////////////
	// Idle Function
	void idle() {}

	//////////////////////////////////////////////////////////////////////////
	// Initial View Setting
	const Vector3f viewEye() const;
	const Vector3f viewRef() const;
	const Vector3f viewUp() const;
	const float viewTheta() const;

	const Vector3f viewEye_Side() const { return Vector3f(1.16162, 0.722146, 2.35805);  }
	const Vector3f viewRef_Side() const { return Vector3f(0.306223, 0.391154, 0.263992); }
	const Vector3f viewUp_Side() const  { return Vector3f(-0.0775126, 0.989164, -0.124687); }
	const float viewTheta_Side() const  { return 45; }

public : 	
};

#endif