#ifndef _VIRTUALSPACE_SCENENODE_H_
#define _VIRTUALSPACE_SCENENODE_H_

#include "../OpenGLState/OpenGLState.h"
#include "../IOHandler/IOHandler.h"
#include "../../Principle.h"

namespace virtualspace {

using virtualspace::OpenGLState;

class Scene;
class PRINCIPLEDLL SceneNode {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	SceneNode(IOHandler *_io=NULL);
	virtual ~SceneNode();

	//////////////////////////////////////////////////////////////////////////
	// glBind
	virtual void glBind() const = 0;

	//////////////////////////////////////////////////////////////////////////
	// Idle Function
	virtual void idle() = 0;

	//////////////////////////////////////////////////////////////////////////
	// Resize Screen
	virtual void resizeScreen(unsigned int _width, unsigned int _height) {}

	//////////////////////////////////////////////////////////////////////////
	// OpenGL State
	void setOpenGLState(const OpenGLState &state)	{ openGLState = state;	}
	const OpenGLState &getOpenGLState() const		{ return openGLState;	}

public : 
	//////////////////////////////////////////////////////////////////////////
	// Get Handler
	IOHandler *getHandler() const { return ioHandler; }

	//////////////////////////////////////////////////////////////////////////
	// Friends
	friend class Scene;

protected : 
	IOHandler *ioHandler;	
	OpenGLState openGLState;
};

};		// namespace virtualspace

#endif

