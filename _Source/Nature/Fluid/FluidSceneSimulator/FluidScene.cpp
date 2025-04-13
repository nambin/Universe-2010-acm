#include <stdafx.h>
#include <Basis/FileIO/FileIO.h>
#include <Basis/Image/Image.h>
#include <Nature/Fluid/FluidScene/FluidScene.h>
#include <Principle/Virtualspace/Scene/Scene.h>

using namespace fluidspace;
using namespace virtualspace;

FluidScene::FluidScene() : integrator(NULL), fluid(NULL), objects(NULL), isBeingSimulated(false) {
	set(1.0f/30.0f);
	ioHandler = new FluidScene_IOHandler(this);
}

FluidScene::FluidScene(real _dt) : integrator(NULL), fluid(NULL), objects(NULL), isBeingSimulated(false) {
	set(_dt);
	ioHandler = new FluidScene_IOHandler(this);
}

FluidScene::~FluidScene() { release(); }

void FluidScene::release() {
	if(integrator)	delete integrator;	integrator = NULL;
	if(fluid)		delete fluid;		fluid = NULL;
	if(objects)		delete objects;		objects = NULL;
}

void FluidScene::set(real _dt) {
	endFrame = 1.e+4;

	frame = 0;
	time = 0;
	dt = _dt;
	isPause = true;
	saveBinaries = saveRenderObjectsOnly = saveImages = false;
}

void FluidScene::idle() {
	if(isPause == false)
		isBeingSimulated = true;

	if(frame >= endFrame) 
		isPause = true;	

	if(!isPause) 
		updateOneStep();

	if(isBeingSimulated == true && isPause == true) {
		post_process();
		isBeingSimulated = false;
	}	
}

void FluidScene::post_process() const {
	integrator->post_process();
	fluid->post_process();
	if(objects) objects->post_process();
}

void FluidScene::saveBinary() const {
	char filename[512];

	sprintf_s( filename, (output_folder + get_concat_filename("FluidBinary/FluidScene", frame)).c_str(), frame );
	std::ofstream ofs(StringLibrary::string2wstring(filename).c_str(), ios::binary);

	::exportBinary(ofs, frame);
	::exportBinary(ofs, endFrame);
	::exportBinary(ofs, dt);
	::exportBinary(ofs, time);

	if(fluid) fluid->saveBinary(frame);
}

void FluidScene::saveRenderObject() const {
	if(fluid) fluid->saveRenderObject(frame);
}

void FluidScene::loadBinary(size_t _frame) {
	char filename[512];

	sprintf_s( filename, (output_folder + get_concat_filename("FluidBinary/FluidScene", _frame)).c_str(), _frame );
	std::ifstream ifs(StringLibrary::string2wstring(filename).c_str(), ios::binary);
	if(!ifs)
		throwError("Error in FluidScene::loadBinary(size_t _frame)");

	::importBinary(ifs, frame);
	::importBinary(ifs, endFrame);
	::importBinary(ifs, dt);
	::importBinary(ifs, time);

	if(fluid) fluid->loadBinary(frame);
}

void FluidScene::saveImage() const {
	char filename[256];
	sprintf_s( filename, (output_folder + get_concat_filename("Fluid/Fluid", frame) + ".tga").c_str(), frame );

	unsigned int width = ioHandler->getScene()->getCamera()->getScreenWidth();
	unsigned int height = ioHandler->getScene()->getCamera()->getScreenHeight();

	glutSwapBuffers();
	Image::saveGLFrameBufferAsTGA(filename, width, height, GL_RGBA);
	glutSwapBuffers();
}

