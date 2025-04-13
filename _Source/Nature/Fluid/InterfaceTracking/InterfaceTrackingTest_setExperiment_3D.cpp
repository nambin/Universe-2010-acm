#include <stdafx.h>
#include <Nature/Fluid/InterfaceTracking/InterfaceTrackingTest.h>
#include <Principle/Geometry/Levelset/TVD_RungeKutta.h>

using namespace geometryspace;
using namespace fluidspace;

void InterfaceTrackingTest::setExperiment_Translation_3D() {
	log_info.set_filename("InterfaceTracking _ Translation 3D");
	test_name = "Translation 3D";

	release();

	unsigned int long_dir = 1;
	grid_index_type long_res = 200, short_res = 30;
	size_t L = 10 , step = 301; 

	real t = 1;
	real dt = t / real(step);

	real h = L / real(long_res);
	real r = short_res * h * (0.5 - 0.2);
	real v = (L - 5.0 * r) / t;		

	Vector3<grid_index_type> dim;		
	Vector3<real> velocity(0);
	Vector3<real> center;

	for(unsigned int i=0;i<3;++i) {
		dim[i] = short_res;
		center[i] = short_res * h * 0.5f;		
	}

	dim[long_dir] = long_res;
	center[long_dir] = r * 2.5;
	velocity[long_dir] = - v;

	//////////////////////////////////////////////////////////////////////////
	// Geometry & Velocity Field	
	SphereDistanceField3D<real> circle(center, r);
	vel = new TranslationVelocityField3D<real>(velocity);

	//////////////////////////////////////////////////////////////////////////
	// Frame
	FluidSimulator::set(dt, 1.0e+6);
	FluidSimulator::setEndFrame(unsigned int(step));

	set_InterfaceTracking(dim[0], dim[1], dim[2], h, circle);
	set_InterfaceIntegrator();	
}

void InterfaceTrackingTest::setExperiment_Zalesak_3D() {
	log_info.set_filename("InterfaceTracking _ Zalesak 3D");
	test_name = "Zalesak 3D";

	release();

	size_t n_rotation = 1;

	size_t res = 64 * 1;
	size_t endTime = 200; 
	real dt = 200.0 / 628.0;

	//res = 100;
	//endTime = 200; 
	//dt = 200.0 / (628.0 * 1);

	// Levelset 
	grid_index_type nX, nY, nZ; nX = nY = nZ = grid_index_type(res);
	real L = 100;
	real h = L / real(res);

	//////////////////////////////////////////////////////////////////////////
	// Geometry & Velocity Field
	ZalesakDistanceField3D<real> zalesak(L);	
	vel = new ZalesakVelocityField3D<real>(endTime, L);

	//////////////////////////////////////////////////////////////////////////
	// Frame
	FluidSimulator::set(dt, 1.0e+6);
	FluidSimulator::setEndFrame(endTime / dt * n_rotation + 1);	

	//set_InterfaceTracking(nX, nY, nZ, h, zalesak);

	GridDistanceField3D<real> bunny; bunny.load("../Data/bunny_phi"); 
	bunny.set_scale(L*0.45);
	bunny.set_rigid_transformation(arithmetic::translate(Vector3<real>(0.5*L,0.75*L,0.5*L)));

	set_InterfaceTracking(nX, nY, nZ, h, bunny);
	set_InterfaceIntegrator();
}

void InterfaceTrackingTest::setExperiment_Enright_3D() {
	log_info.set_filename("InterfaceTracking _ Enright Test 3D");
	test_name = "Enright Test 3D";

	release();

	grid_index_type res = 100;
	unsigned int frame = 600;

	grid_index_type nX, nY, nZ;
	nX = nY = nZ = res;

	real L = 1.0f;
	real h = L / nX;	

	SphereDistanceField3D<real> sphere(Vector3<real>(0.35*L), 0.15*L);
	vel = new EnrightVelocityField3D<real>(L);

	set(3.0f / float(frame), 1.0e+6);
	setEndFrame(frame+1);	

	set_InterfaceTracking(nX, nY, nZ, h, sphere);
	set_InterfaceIntegrator();

	tracking->eye	= Vector3f(0.5f * L, - L * 1.25f, 0.5f * L);
	tracking->ref	= Vector3f(0.5f * L,		  0	, 0.5f * L);
	tracking->up	= Vector3f(0,0,1);
	tracking->theta = 45;
}

void InterfaceTrackingTest::setExperiment_Object_3D(const std::string& obj_file) {
	log_info.set_filename("InterfaceTracking _ Object 3D");
	test_name = "Object 3D";

	release();

	size_t n_rotation = 1;

	size_t res = 64 * 2;
	size_t endTime = 200; 
	real dt = 200.0 / 628.0;

	//res = 100;
	//endTime = 200; 
	//dt = 200.0 / (628.0 * 1);

	// Levelset 
	grid_index_type nX, nY, nZ; nX = nY = nZ = grid_index_type(res);
	real L = 2;
	real h = L / real(res);

	//////////////////////////////////////////////////////////////////////////
	// Geometry & Velocity Field
	ObjectDistanceField3D<real> field(obj_file.c_str()); field.set_rigid_transformation(arithmetic::translate(Vector3<real>(L*0.5,L*0.5,L*0.5)));
	vel = new ZalesakVelocityField3D<real>(endTime, L);

	//////////////////////////////////////////////////////////////////////////
	// Frame
	FluidSimulator::set(dt, 1.0e+6);
	FluidSimulator::setEndFrame(endTime / dt * n_rotation + 1);	

	set_InterfaceTracking(nX, nY, nZ, h, field);
	set_InterfaceIntegrator();
}