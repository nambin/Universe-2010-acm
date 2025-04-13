#include <stdafx.h>
#include <Principle/VirtualSpace/Scene/Scene.h>

using namespace virtualspace;

Scene::Scene() {
	init();
}

Scene::~Scene() {
	release();
}

void Scene::init() {
	
	ioHandler = new Scene_IOHandler(this);
	clearColor = RGBA::black();

	camera = new Camera();
	camera->resizeScreen(Scene::initWindowWidth(), Scene::initWindowHeight());

	lightSize = 0.1f;
	idxSelectedObject = -1;	

	useTBB = false;
}

void Scene::release() {
	if(ioHandler)
		delete ioHandler; 
	ioHandler = NULL;
}

void Scene::resizeScreen(unsigned int _width, unsigned int _height) { 
	camera->resizeScreen(_width, _height);  
	for(size_t i=0;i<sceneObjects.size();++i)
		sceneObjects[i]->resizeScreen(_width, _height);
}

void Scene::setLightSize(float size) {
	lightSize = size;
}

void Scene::push_SceneLight(SceneLight *light) {
	if(light == NULL)
		return;

	light->getHandler()->setScene(this);
	sceneLights.push_back(light);
}

void Scene::push_SceneLight(const Light &light) {
	SceneLight *scenelight = new SceneLight(light);
	push_SceneLight(scenelight);
}

void Scene::push_SceneObject(SceneObject *object) {
	if(object != NULL) {
		object->getHandler()->setScene(this);
		object->resizeScreen(camera->getScreenWidth(), camera->getScreenHeight());

		sceneObjects.push_back(object);		
	}
	idxSelectedObject = 0;
}

void Scene::set_SceneLightVisible(bool isvisible) {
	for(size_t i=0;i<sceneLights.size();++i)
		sceneLights[i]->setVisible(isvisible);
}

void Scene::selectNextSceneObject() const {
	++idxSelectedObject;
	if(idxSelectedObject >= (int)sceneObjects.size()) 
		idxSelectedObject = 0;
}

#ifdef USE_TBB
void Scene::enableUseTBB(unsigned int nThreads) { 
	useTBB = true; 
	tbbprogram.initialze(nThreads);
}
#endif