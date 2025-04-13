#include <stdafx.h>
#include <Nature/Fluid/FluidScene/FluidScene.h>

using namespace std;
using namespace fluidspace;

void FluidScene::updateOneStep() {
	if(saveBinaries) {
		if(saveRenderObjectsOnly && frame <= saving_frame_renderObject)	saveRenderObject();
		else if(frame <= saving_frame)									saveBinary();
	}
	if(saveImages) {
		saveImage();
	}

	update(dt);	
	glutPostRedisplay();
}

void FluidScene::update(real dt) {
	if(fluid)	
		integrator->update(dt, objects, fluid);	

	time += dt;
	frame++;
}

void FluidScene::initFluidSimulator(FluidSimulator* sim) {
	if(fluid) delete fluid; 
	fluid = sim; 

	if(fluid) {
		fluid->resizeScreen(getHandler()->getScene()->getCamera()->getScreenWidth(), getHandler()->getScene()->getCamera()->getScreenHeight());
		fluid->saveImages = saveImages;
		fluid->saveBinaries = saveBinaries;
	}	
}

void FluidScene::initFluidScene_Objects(FluidScene_Objects* objs) {
	if(objects) delete objects; 
	objects = objs;
}

void FluidScene::initFluidSceneIntegrator(FluidSceneIntegrator* itg) {
	if(integrator) delete integrator; 
	integrator = itg;
}