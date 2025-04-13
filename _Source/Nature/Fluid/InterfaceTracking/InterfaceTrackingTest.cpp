#include <stdafx.h>
#include <Nature/Fluid/InterfaceTracking/InterfaceTrackingTest.h>
#include <Principle/Geometry/Levelset/TVD_RungeKutta.h>

using namespace geometryspace;
using namespace fluidspace;

InterfaceTrackingTest::InterfaceTrackingTest() {	
	vel = NULL;
	tracking = NULL;
	integrator = NULL;
	ioHandler = new InterfaceTrackingTest_IOHandler(this);

	dBoundingBox = dLiquid = dFilledMesh = true;
	dGrid = dPhi = dAuxiliary = false;	
}

void InterfaceTrackingTest::release() {
	if(vel) delete vel; vel = NULL;
	if(tracking) delete tracking; tracking = NULL;
	if(integrator) delete integrator; integrator = NULL;
}

void InterfaceTrackingTest::update(real _dt) {
	if(vel->getDimension() != tracking->getDimension()) {
		log_info << "Error in InterfaceTrackingTest::update(real _dt)" << endl;
		exit(0);
	}

	log_info << endl << "Frame " << frame << " at t = " << time << endl;

	Timer timer;

	integrator->update(*tracking, *vel, _dt);		timer.logTime("Interface Integration", log_info);	timer.reset();
	vel->addTime(_dt);

	if(saveImages || tracking->getDimension() == 2) {
		tracking->createDrawableInterface();			
		timer.logTime("Create Drawable Interface", log_info);
	}	

	frame++; 
	time += _dt;

	char str[256];
	sprintf_s(str, "Liquid Volume : %30.20g", tracking->getTotalVolume());
	log_info << str << endl;
}

void InterfaceTrackingTest::saveBinary(std::size_t frame) const {
	char filename[512];
	sprintf_s( filename, (output_folder + get_concat_filename("FluidBinary/%s _ " + tracking->name_class(), frame)).c_str(), test_name.c_str(), frame );

	tracking->save_compressed(filename);
}

void InterfaceTrackingTest::post_process() const {
	if(tracking)	tracking->post_process();
	if(integrator)	integrator->post_process();
}