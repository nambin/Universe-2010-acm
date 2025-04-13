#include <stdafx.h>
#include <Nature/MassSpring/ClothScene/ClothScene.h>
#include <Basis/Image/Image.h>
#include <Principle/Virtualspace/Scene/Scene.h>

using namespace springspace;
using namespace virtualspace;

ClothScene::ClothScene() : mesh(NULL), objects(NULL), integrator(NULL) {
	set(1.0f/30.0f);
	ioHandler = new ClothScene_IOHandler(this);

	dSurface = true;
	dParticle = dStretchSpring = dShearSpring = dBendingSpring = false;
	dObjects = true;
}

ClothScene::ClothScene(real _dt) : mesh(NULL), objects(NULL), integrator(NULL) {
	set(_dt);
	ioHandler = new ClothScene_IOHandler(this);

	dSurface = true;
	dParticle = dStretchSpring = dShearSpring = dBendingSpring = false;
	dObjects = true;
}

ClothScene::~ClothScene() { release(); }

void ClothScene::release() {
	if(mesh)		delete mesh;		mesh = NULL;
	if(objects)		delete objects;		objects = NULL;
	if(integrator)	delete integrator;	integrator = NULL;
}

void ClothScene::set(real _dt) {
	endFrame = 1.e+4;

	frame = 0;
	time = 0;
	dt = _dt;
	isPause = true;
	saveImages = false;
}

void ClothScene::idle() {
	if(frame >= endFrame)
		isPause = true;

	if(!isPause) 
		updateOneStep();		
}

void ClothScene::saveImage() const {
	char filename[256];
	sprintf_s( filename, (output_folder + "Cloth/Cloth _ %4d" + ".tga").c_str(), frame );	

	unsigned int width = ioHandler->getScene()->getCamera()->getScreenWidth();
	unsigned int height = ioHandler->getScene()->getCamera()->getScreenHeight();

	Image::saveGLFrameBufferAsTGA(filename, width, height, GL_RGBA);
}

