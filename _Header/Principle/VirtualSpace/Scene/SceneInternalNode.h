#ifndef _VIRTUALSPACE_SCENE_INTERNALNODE_H_
#define _VIRTUALSPACE_SCENE_INTERNALNODE_H_

#include "SceneNode.h"
#include "SceneStateNode.h"
#include "ScenePrimitive.h"

namespace virtualspace {

class Scene;
class PRINCIPLEDLL SceneInternalNode : public SceneNode {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	SceneInternalNode();
	~SceneInternalNode();

	//////////////////////////////////////////////////////////////////////////
	// glBind
	void glBind() const;

	//////////////////////////////////////////////////////////////////////////
	// Idle Function
	void idle();	

private : 
	std::vector<SceneStateNode*>	nodeStates;
	std::vector<ScenePrimitive*>	primitives;
	std::vector<SceneInternalNode*>	children;
};

};	// namespace virtualspace

#endif

