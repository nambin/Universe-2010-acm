#include <stdafx.h>
#include <Principle/VirtualSpace/Scene/SceneTransform.h>

using namespace virtualspace;

void SceneTransform::glBind() const {
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glTransform();
}

void SceneTransform::glUnbind() const {
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}

void SceneTransform::idle() {
}