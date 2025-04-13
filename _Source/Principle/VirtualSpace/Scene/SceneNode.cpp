#include <stdafx.h>
#include <Principle/VirtualSpace/Scene/SceneNode.h>

using namespace virtualspace;

SceneNode::SceneNode(IOHandler *_io) {
	ioHandler = _io;	
}

SceneNode::~SceneNode() {
	if(ioHandler)
		delete ioHandler;
}

void SceneNode::glBind() const {	
	openGLState.glBind();
}
