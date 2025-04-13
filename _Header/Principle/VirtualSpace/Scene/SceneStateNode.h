#ifndef _VIRTUALSPACE_SCENE_STATENODE_H_
#define _VIRTUALSPACE_SCENE_STATENODE_H_

#include "SceneNode.h"

namespace virtualspace {

class Scene;
class PRINCIPLEDLL SceneStateNode : public SceneNode {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	SceneStateNode(IOHandler *_io=NULL);
	virtual ~SceneStateNode();

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	virtual void glBind() const = 0;
	virtual void glUnbind() const = 0;

	//////////////////////////////////////////////////////////////////////////
	// Idle Function
	virtual void idle() = 0;
};

};	// namespace virtualspace

#endif