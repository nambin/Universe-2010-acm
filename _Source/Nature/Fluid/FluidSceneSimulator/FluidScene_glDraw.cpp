#include <stdafx.h>
#include <Nature/Fluid/FluidScene/FluidScene.h>

using namespace std;
using namespace fluidspace;

void FluidScene::glBind() const {
	SceneObject::glBind();

	OpenGLState state = openGLState;
	state.glBind();

	glDraw();	
}

void FluidScene::glDraw() const {		
	if(fluid)	fluid->glBind();
	if(objects)	objects->glDraw();
}

void FluidScene::resizeScreen(unsigned int width, unsigned int height) {
	if(fluid) fluid->resizeScreen(width, height);
}

//////////////////////////////////////////////////////////////////////////
// Initial View Setting
const Vector3f FluidScene::viewEye() const { 
	if(fluid)	return fluid->viewEye();
	return Vector3f(0,0,1);
}

const Vector3f FluidScene::viewRef() const { 
	if(fluid)	return fluid->viewRef();
	return Vector3f(0,0,0);
}

const Vector3f FluidScene::viewUp() const  { 
	if(fluid)	return fluid->viewUp();
	return Vector3f(0,1,0);
}

const float FluidScene::viewTheta() const  { 
	if(fluid)	return fluid->viewTheta();
	return 45;
}
