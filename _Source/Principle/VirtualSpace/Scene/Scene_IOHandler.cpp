#include <stdafx.h>
#include <Principle/VirtualSpace/Scene/Scene.h>

using namespace virtualspace;

Scene_IOHandler::Scene_IOHandler(Scene *s)  {
	setScene(s);
}

Scene_IOHandler::~Scene_IOHandler() {
}

void Scene_IOHandler::keyboardHandler(unsigned char key, int x, int y) {
	scene->camera->getHandler()->keyboardHandler(key, x, y);

	//if(key == ' ') {
	//	scene->selectNextSceneObject();
	//}
	//else {
		scene->sceneObjects[scene->idxSelectedObject]->getHandler()->keyboardHandler(key, x, y);
	//}
}

void Scene_IOHandler::mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y) {	
	scene->camera->getHandler()->mouseHandler(button, state, x, y);

	if(scene->idxSelectedObject >= 0 && scene->idxSelectedObject < static_cast<int>(scene->sceneObjects.size())) {
		scene->sceneObjects[scene->idxSelectedObject]->getHandler()->mouseHandler(button, state, x, y);
	}
}

void Scene_IOHandler::mouseMotionHandler(int x, int y) {
	scene->camera->getHandler()->mouseMotionHandler(x, y);
	scene->sceneObjects[scene->idxSelectedObject]->getHandler()->mouseMotionHandler(x,y);
}

void Scene::idle() {
	if(idxSelectedObject > -1 && idxSelectedObject < (int)sceneObjects.size()) 
		sceneObjects[idxSelectedObject]->idle();
}