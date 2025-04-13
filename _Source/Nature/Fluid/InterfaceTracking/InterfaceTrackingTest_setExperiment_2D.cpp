#include <stdafx.h>
#include <Nature/Fluid/InterfaceTracking/InterfaceTrackingTest.h>
#include <Principle/Geometry/Levelset/TVD_RungeKutta.h>

using namespace geometryspace;
using namespace fluidspace;

void InterfaceTrackingTest::setExperiment_Translation_2D() {
	log_info.set_filename("InterfaceTracking _ Translation 2D");
	test_name = "Translation 2D";

	release();

	grid_index_type nX = 200 , nY = 15;
	size_t L = 1;
	unsigned int step = 801; 
	bool swap_dir = false;
	
	real t = 1;
	real dt = t / real(step);
	
	real h = L / real(nX);
	real r = nY * h * (0.5 - 0.2);
	real v = (L - 5.0 * r) / t;		

	Vector2<real> velocity(v,0);
	Vector2<real> center(r * 2.5, nY * h * 0.5);
	if(swap_dir) {
		swap(nX, nY); swap(velocity[0], velocity[1]); swap(center[0], center[1]);
	}

	//////////////////////////////////////////////////////////////////////////
	// Geometry & Velocity Field	
	CircleDistanceField2D<real> circle(center, r);
	vel = new TranslationVelocityField2D<real>(velocity);

	//////////////////////////////////////////////////////////////////////////
	// Frame
	FluidSimulator::set(dt, 1.0e+6);
	FluidSimulator::setEndFrame(step);

	set_InterfaceTracking(nX, nY, h, circle);
	set_InterfaceIntegrator();
}

void InterfaceTrackingTest::setExperiment_Zalesak_2D() {
	log_info.set_filename("InterfaceTracking _ Zalesak 2D");
	test_name = "Zalesak 2D";

	release();

	size_t n_rotation = 1;

	size_t res = 128;
	size_t endTime = 200; 
	real dt = 200.0 / 628.0;

	//size_t res = 200;
	//size_t endTime = 200; 
	//real dt = 200.0 / (628.0 * 2);

	// Levelset 
	grid_index_type nX, nY; nX = nY = grid_index_type(res);	//nX = res * 1.5;
	real L = 100;
	real h = L / real(res);

	//////////////////////////////////////////////////////////////////////////
	// Geometry & Velocity Field
	ZalesakDistanceField2D<real> zalesak(L);	
	vel = new ZalesakVelocityField2D<real>(endTime, L);

	//////////////////////////////////////////////////////////////////////////
	// Frame
	FluidSimulator::set(dt, 1.0e+6);
	FluidSimulator::setEndFrame(endTime / dt * n_rotation + 0);
	//FluidSimulator::setEndFrame(30);

	set_InterfaceTracking(nX, nY, h, zalesak);
	set_InterfaceIntegrator();
}

