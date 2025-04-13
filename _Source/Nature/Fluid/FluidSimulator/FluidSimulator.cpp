#include <stdafx.h>
#include <Basis/Image/Image.h>
#include <Basis/FileIO/FileIO.h>
#include <Nature/Fluid/FluidSimulator/FluidSimulator.h>
#include <Principle/VirtualSpace/Scene/Scene.h>

using namespace fluidspace;
using namespace virtualspace;
using namespace basis;

FluidSimulator::FluidSimulator() : objects(NULL) {
	set(0,0);
	image_filename = "Fluid";
}

FluidSimulator::FluidSimulator(real _dt, real _cfl) : objects(NULL) {
	set(_dt, _cfl);
	image_filename = "Fluid";
}

FluidSimulator::~FluidSimulator() {}

void FluidSimulator::set(real _dt, real _cfl) {
	endFrame = 1.0e+4;

	frame = 0;
	time = 0;
	dt = dt0 = _dt;
	maxCFL = _cfl;
	isPause = true;
	saveBinaries = saveImages = false;
}

void FluidSimulator::glBind() const {
	SceneObject::glBind();
}

void FluidSimulator::idle() {
	if(frame >= endFrame)
		isPause = true;

	if(!isPause) {
		updateOneStep();		
		glutPostRedisplay();
	}
}

void FluidSimulator::updateOneStep() {
	if(saveBinaries)
		saveBinary(frame);

	if(saveImages)
		saveImage();

	update(dt0);		
}

void FluidSimulator::loadBinary(size_t frame) {
	basis::throwError("Not Implemented _ FluidSimulator::loadinary()");
}

void FluidSimulator::saveBinary(size_t frame) const {
	basis::throwError("Not Implemented _ FluidSimulator::saveBinary()");
}

void FluidSimulator::saveRenderObject(size_t frame) const {
	basis::throwError("Not Implemented _ FluidSimulator::saveRenderObject()");
}

void FluidSimulator::importBinary(std::ifstream& ifs) {
	::importBinary(ifs, frame);
	::importBinary(ifs, endFrame);
	::importBinary(ifs, dt);
	::importBinary(ifs, dt0);
	::importBinary(ifs, time);
	::importBinary(ifs, maxCFL);
}

void FluidSimulator::exportBinary(std::ofstream& ofs) const {
	::exportBinary(ofs, frame);
	::exportBinary(ofs, endFrame);
	::exportBinary(ofs, dt);
	::exportBinary(ofs, dt0);
	::exportBinary(ofs, time);
	::exportBinary(ofs, maxCFL);
}

void FluidSimulator::saveImage() const {
	char filename[256];
	sprintf_s( filename, (output_folder + get_concat_filename("Fluid/" + image_filename, frame) + ".tga").c_str(), frame );

	unsigned int width = ioHandler->getScene()->getCamera()->getScreenWidth();
	unsigned int height = ioHandler->getScene()->getCamera()->getScreenHeight();

	glutSwapBuffers();
	Image::saveGLFrameBufferAsTGA(filename, width, height, GL_RGBA);
	glutSwapBuffers();
}



