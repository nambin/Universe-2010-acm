#ifndef _VIRTUALSPACE_SCENE_PRIMITIVE_H_
#define _VIRTUALSPACE_SCENE_PRIMITIVE_H_

#include "SceneNode.h"

namespace virtualspace {

class Scene;
class PRINCIPLEDLL ScenePrimitive : public SceneNode {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	ScenePrimitive(IOHandler *_io=NULL);
	virtual ~ScenePrimitive();

	//////////////////////////////////////////////////////////////////////////
	// glBind
	virtual void glBind() const = 0;

	//////////////////////////////////////////////////////////////////////////
	// Idle Function
	virtual void idle() = 0;
};

};		// namespace virtualspace

#endif