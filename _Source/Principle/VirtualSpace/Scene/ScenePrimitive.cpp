#include <stdafx.h>
#include <Principle/VirtualSpace/Scene/ScenePrimitive.h>

using namespace virtualspace;

ScenePrimitive::ScenePrimitive(IOHandler *_io) : SceneNode(_io) {}

ScenePrimitive::~ScenePrimitive() {}

void ScenePrimitive::glBind() const {
	SceneNode::glBind();
}