#include <stdafx.h>
#include <Principle/VirtualSpace/RenderingPipeline/RenderingPipelineSetting.h>
#include <Principle/VirtualSpace/DrawLibrary/glDrawLibrary.h>
#include <Principle/VirtualSpace/Scene/Scene.h>

using namespace virtualspace;

void Scene::glSetup() const {
	glClearColor(clearColor[0], clearColor[1], clearColor[2], clearColor[3]);	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void Scene::glDraw() const {

	glSetup();

	// Camera Setup
	camera->glBind();

	// Lights Setup
	Light::initLights();
	for(unsigned int i=0;i<sceneLights.size();++i) {
		sceneLights[i]->glBind();	
		sceneLights[i]->glDraw(lightSize);
	}	

	// Rendering
	renderScene();

	// Draw View Orientation	
	GlDrawLibrary::glDrawViewOrientation(camera->getScreenWidth(), camera->getScreenHeight(), 0.10f, camera->getModelViewMatrix().getptr());	
}

void Scene::renderScene() const {
	for(unsigned int i=0;i<sceneObjects.size();++i) 
		sceneObjects[i]->glBind();
}
