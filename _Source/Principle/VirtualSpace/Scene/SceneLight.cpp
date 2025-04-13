#include <stdafx.h>
#include <Principle/VirtualSpace/Scene/SceneLight.h>

using namespace virtualspace;

SceneLight::SceneLight() : SceneStateNode(), visible(false) {
	ioHandler = new SceneLight_IOHandler(this);
}

SceneLight::SceneLight(const Light &light) : SceneStateNode(), visible(false) {
	setLight(light);
	ioHandler = new SceneLight_IOHandler(this);
}

SceneLight::~SceneLight() {}

void SceneLight::glBind() const {
	light.glBindLight();
}

void SceneLight::glUnbind() const {
	light.glUnbindLight();
}

void SceneLight::glDraw(float size) const {
	if(!visible)
		return;
	light.glDraw(size);
}

void SceneLight::idle() {
}