#ifndef _VIRTUALSPACE_SCENE_H_
#define _VIRTUALSPACE_SCENE_H_

#include <vector>
#include "../../../Basis/TBBProgramming/TBBProgramming.h"
#include "../../Principle.h"
#include "../IOHandler/IOHandler.h"
#include "../Camera/Camera.h"
#include "../Scene/SceneLight.h"
#include "../Scene/SceneObject.h"

namespace virtualspace {

using namespace std;

class PRINCIPLEDLL Scene_IOHandler : public IOHandler {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Scene_IOHandler(Scene *s);
	~Scene_IOHandler();

	//////////////////////////////////////////////////////////////////////////
	// Handler
	void keyboardHandler(unsigned char key, int	x, int y);
	void mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y);
	void mouseMotionHandler(int x, int y);
};

class PRINCIPLEDLL Scene {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Scene();
	~Scene();

	void init();
	void release();

	void push_SceneLight(SceneLight *);
	void push_SceneLight(const Light &);
	void set_SceneLightVisible(bool isVisible = true);
	void setLightSize(float);

	void push_SceneObject(SceneObject *);
	void selectNextSceneObject() const;

	//////////////////////////////////////////////////////////////////////////
	// Handler
	void keyboardHandler(unsigned char k, int x, int y) { ioHandler->keyboardHandler(k, x, y); glutPostRedisplay(); }
	void mouseHandler(IOHandler::IO_MOUSE_BUTTON b, IOHandler::IO_MOUSE_STATE s, int x, int y) { ioHandler->mouseHandler(b, s, x, y); glutPostRedisplay(); }
	void mouseMotionHandler(int x, int y) { ioHandler->mouseMotionHandler(x,y); glutPostRedisplay(); }

	void resizeScreen(unsigned int _width, unsigned int _height);

	void idle();

	//////////////////////////////////////////////////////////////////////////
	// Getter
	Camera *getCamera() const { return camera; }

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	void glSetup() const;
	void glDraw() const;
	void renderScene() const;

#ifdef USE_TBB
	//////////////////////////////////////////////////////////////////////////
	// TBB
	virtual void enableUseTBB(unsigned int nThreads);
#endif

	//////////////////////////////////////////////////////////////////////////
	// Friends 
	friend Scene_IOHandler;

public : 
	static unsigned int initWindowWidth()	{ return 800;	}
	static unsigned int initWindowHeight()	{ return 800;	}
	static unsigned int InitWindowX()		{ return 100;	}
	static unsigned int InitWindowY()		{ return 100;	}

private : 
	Scene_IOHandler *			ioHandler;

	RGBA						clearColor;

	mutable Camera *			camera;
	
	float						lightSize;
	std::vector<SceneLight*>	sceneLights;

	std::vector<SceneObject*>	sceneObjects;		
	mutable int					idxSelectedObject;

	//////////////////////////////////////////////////////////////////////////
	// Variable for TBB
	bool					useTBB;

#ifdef USE_TBB
	basis::TBBProgramming	tbbprogram;
#endif
};

};	// namespace virtualspace

#endif