#include <stdafx.h>
#include "SRL_Experiment.h"

using namespace std;
using namespace geometryspace;
using namespace rigidspace;
using namespace fluidspace;

typedef fluidspace::real T;

void SRL_FluidSimulation_3D_Demo_05(Scene* scene) {
#ifdef USE_OPENMP	
	const int maxthreads = omp_get_max_threads();
	const int numThreads = maxthreads;

	omp_set_num_threads(numThreads);
	log_info << endl << "OpenMP " << numThreads << " Threads are initiated" << endl << endl;
#endif

	//const real g = - 9.81 * 1.0f, sigma = 0.000;

	//const int ratio = pow(2.0, 4);
	//const int nX = 12 * ratio, nY = 9 * ratio, nZ = 12 * ratio;
	//const real h = 1.0 / double(nY);	
	//const real Lx = nX*h, Ly = nY*h, Lz = nZ*h;

	//const real water_center_y = Ly*0.85, water_radius = Ly*0.06;
	//const real water_vel_y = Ly * 500.0f;

	//real cube_h = Ly * 0.4;
	//real cube_width = Ly * 0.2;

	//const int endFrame = 99;
	//const real dt = 1.0 / 40000.0;	

	//const real g = - 9.81 * 0.0f, sigma = 0.0;

	//const int ratio = 96;
	//const int nX = 2 * ratio, nY = 1 * ratio, nZ = 2 * ratio;
	//const real h = 1.0 / double(nY);	
	//const real Lx = nX*h, Ly = nY*h, Lz = nZ*h;

	//const real water_center_x = Lx*0.49, water_center_z = Lz*0.52;

	//const real water_center_y = Ly*0.85, water_radius = Ly * 0.09;
	//const real water_vel_y = Ly * 300.0f;

	//real cube_h = Ly * 0.4;
	//real cube_width = Ly * 0.3;

	//const int endFrame = 99;
	//const real dt = 1.0 / 20000.0;	
	
	const real g = - 9.81 * 0.0f, sigma = 0.0;

	const int ratio = 96;
	const int nX = 2 * ratio, nY = 1 * ratio, nZ = 2 * ratio;
	const real h = 1.0 / double(nY);	
	const real Lx = nX*h, Ly = nY*h, Lz = nZ*h;

	const real water_center_x = Lx*0.49, water_center_z = Lz*0.52;

	const real water_center_y = Ly*0.85, water_radius = Ly * 0.09;
	const real water_vel_y = Ly * 300.0f;

	real cube_h = Ly * 0.4;
	real cube_width = Ly * 0.3;

	const int endFrame = 99;
	const real dt = 1.0 / 20000.0;		

	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////
	const real cfl = 3.5, cfl_interface = 1.5;

	const real narrowband_srl = 5;
	const unsigned int redist_freq = 15;
	const unsigned int redist_min = 5, redist_max = 50;	

	FluidScene *fluidscene = new FluidScene(dt);
	fluidscene->initFluidSceneIntegrator(new FluidSceneIntegrator_StaticRigidObject(dt));	
	fluidscene->setEndFrame(endFrame);
	fluidscene->enableSaveBinary(endFrame);
	fluidscene->enableSaveRenderObjectsOnly(endFrame);
//	fluidscene->enableSaveImage();		
	scene->push_SceneObject(fluidscene);

	Liquid3D *liquid3D = new Liquid3D(nX, nY, nZ, h, dt, cfl, cfl_interface);
	liquid3D->release();

	////////////////////////////////////////////////////////////////////////	
	// Interface Tracking
	liquid3D->tracking = new InterfaceTracking_SRL3D(nX, nY, nZ, h, narrowband_srl);

	////////////////////////////////////////////////////////////////////////
	// CFD Integrator
	CFDIntegrator_Liquid* itg = new CFDIntegrator_Liquid();	

	itg->set_CFDAddForce(new CFDAddForce_Gravity(g));	
	itg->set_CFDAdvection(new CFDAdvection_SemiLagrangian<>());

	CFDProjection_Liquid* project = new CFDProjection_FreeSurface_Gibou(false, sigma);

	project->set_InterfaceTracking(liquid3D->tracking);
	project->init_PoissonSolver_GMG(liquid3D->solver, 12); project->set_iteration_number_GMG(3,5);
	project->set_divergence_threshold(1.0e-3);

	itg->set_CFDProjection(project);

	//itg->set_CFDTimeStepRestriction(new CFDTimeStepRestriction(cfl));	
	itg->set_CFDTimeStepRestriction(new CFDTimeStepRestriction_SurfaceTension(cfl, g, sigma));	

	itg->set_CFDVelocityBoundaryCondition(new GridBoundaryCondition_U_FreeSlip(), new GridBoundaryCondition_V_FreeSlip(), new GridBoundaryCondition_W_FreeSlip());

	////////////////////////////////////////////////////////////////////////	
	// Interface Integrator
	InterfaceIntegrator_SRL* itg_srl = new InterfaceIntegrator_SRL();

	SRL_Redistancing* redist = new SRL_Redistancing_PDE_WENO(narrowband_srl, 0.5, true, 0);
	redist->set_redistance_cfl_frequency(redist_freq);	
	redist->set_redistance_min_max_updates(redist_min, redist_max);	

	itg_srl->set_Advection(new CFDAdvection_SemiLagrangian_SRL<CFDTracking_RK2>(true));		
	itg_srl->set_Redistancing(redist);		
	itg_srl->set_TimeStepRestriction(new CFDTimeStepRestriction(cfl_interface));

	Levelset_Extrapolation<real>* extrapolation = new Levelset_Extrapolation_PDE<real>(narrowband_srl, 0.5f);

	LiquidSource* source = NULL;

	itg->set(liquid3D->tracking, itg_srl, extrapolation, source);	
	liquid3D->integrator = itg;

	////////////////////////////////////////////////////////////////////////	
	// Set Phi
	JointDistanceField3D<real> jointfield;
	
	WaterDropDistanceField3D<real> water(Vector3<real>(water_center_x, water_center_y, water_center_z), water_radius, -3.00*Ly);
	jointfield.add_distancefields(&water);
	liquid3D->setPhi(jointfield);

	LiquidSource_Discrete* liquid_source = new LiquidSource_Discrete(NULL);
	liquid_source->setVel(liquid3D->solver, Vector3<real>(0, - water_vel_y, 0), water, 5.0f);

	fluidscene->initFluidSimulator(liquid3D);	

	////////////////////////////////////////////////////////////////////////	
	// Set Object
	FluidScene_Objects3D* object3D = new FluidScene_Objects3D();

	Vector3<real> cube_min(Lx*0.5 - cube_width*0.5, 0, Lz*0.5 - cube_width*0.5);
	Vector3<real> cube_max(Lx*0.5 + cube_width*0.5, cube_h, Lz*0.5 + cube_width*0.5);

	// Cube
	log_info << "Cube " << cube_min << " _ " << cube_max << std::endl;	
	object3D->push_DynamicObject(StaticRigidBody3D<T>(Cube<T>(cube_min, cube_max)));

	//// Sphere
	//log_info << "Sphere " << Vector3<real>(Lx*0.5, 0, Lz*0.5) << " _ " << cube_h << std::endl;	
	//object3D->push_DynamicObject(StaticRigidBody3D<T>(Sphere<T>(Vector3<real>(Lx*0.5, 0, Lz*0.5), cube_h*1.0)));

	fluidscene->initFluidScene_Objects(object3D);
}