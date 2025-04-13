#include <stdafx.h>
#include <Principle/VirtualSpace/Scene/SceneInternalNode.h>

using namespace virtualspace;

SceneInternalNode::SceneInternalNode() : SceneNode(NULL) {}

SceneInternalNode::~SceneInternalNode() {}

void SceneInternalNode::glBind() const {
	SceneNode::glBind();

	for(unsigned int i=0;i<nodeStates.size();++i)
		nodeStates[i]->glBind();

	for(unsigned int i=0;i<primitives.size();++i)
		primitives[i]->glBind();

	for(unsigned int i=0;i<children.size();++i)
		children[i]->glBind();

	for(unsigned int i=0;i<nodeStates.size();++i)
		nodeStates[i]->glUnbind();	
}

void SceneInternalNode::idle() {
}