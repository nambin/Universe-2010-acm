#include <stdafx.h>
#include <Nature/MassSpring/ClothScene/ClothScene.h>

using namespace springspace;

void ClothScene::initDynamicMesh(DynamicMesh* _mesh) {
	if(mesh)
		delete mesh;
	mesh = _mesh;
}

void ClothScene::initClothScene_Objects(ClothScene_Objects* scene) {
	if(objects) 
		delete objects;
	objects = scene;
}

void ClothScene::initClothSceneIntegrator(ClothSceneIntegrator* _itg) {
	if(integrator)
		delete integrator;
	integrator = _itg;
}

void ClothScene::updateOneStep() {
	if(saveImages)
		saveImage();

	update(dt);	
}

void ClothScene::update(real dt) {
	integrator->update(dt, objects, mesh);	
	mesh->postProcessing();
	time += dt;
	frame++;

	glutPostRedisplay();
}
