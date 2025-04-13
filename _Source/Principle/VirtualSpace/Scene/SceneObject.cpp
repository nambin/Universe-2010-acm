#include <stdafx.h>
#include <Principle/VirtualSpace/Scene/SceneObject.h>

using namespace virtualspace;

SceneObject::SceneObject(IOHandler *_io) : ScenePrimitive(_io) {}
SceneObject::~SceneObject() {}

void SceneObject::glBind() const {
	ScenePrimitive::glBind();
	material.glBindMaterial();
}
