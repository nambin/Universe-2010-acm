#include <stdafx.h>
#include "SRL_Experiment.h"

using namespace std;
using namespace geometryspace;
using namespace rigidspace;
using namespace fluidspace;

typedef fluidspace::real T;

void SRL_FluidSimulation_DamBreaking(Scene* scene) {
	const real g = - 9.81, sigma = 0.001;

	const int endFrame = 200;
	const real dt = 1.0 / 30.0;	
	const real cfl = 1.0f, cfl_interface = 1.0f;

	const real narrowband_srl = 4;
	const unsigned int redist_freq = 7;
	const unsigned int redist_min = 4, redist_max = 25;

	const int ratio = 1;
	const int nX = 64*ratio, nY = 96*ratio, nZ = nX; 
	const real h = 1.0 / double(nX);	

	const real Lx = nX*h, Ly = nY*h, Lz = nZ*h;

	FluidScene *fluidscene = new FluidScene(dt);
	fluidscene->initFluidSceneIntegrator(new FluidSceneIntegrator_StaticRigidObject(dt));	
	fluidscene->setEndFrame(endFrame);
	fluidscene->enableSaveBinary(500);
	fluidscene->enableSaveRenderObjectsOnly(500);
	fluidscene->enableSaveImage();		
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

	CFDProjection_Liquid* project = new CFDProjection_FreeSurface_Gibou(true, sigma);
	project->set_InterfaceTracking(liquid3D->tracking);
	project->init_PoissonSolver_GMG(liquid3D->solver, 12); project->set_iteration_number_GMG(3,5);
	project->set_divergence_threshold(1.0e-4);
	itg->set_CFDProjection(project);

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

	Levelset_Extrapolation<real>* extrapolation = NULL;
	extrapolation = new Levelset_Extrapolation_PDE<real>(narrowband_srl, 0.5f);

	itg->set(liquid3D->tracking, itg_srl, extrapolation);	
	liquid3D->integrator = itg;

	////////////////////////////////////////////////////////////////////////	
	// Set Phi
	liquid3D->setPhi(CubeDistanceField3D<real>(Vector3<real>(-10.0f * Lx), Vector3<real>(0.4f*Lx, 0.45*Ly, 0.4f*Lz)));
	fluidscene->initFluidSimulator(liquid3D);	
}

void SRL_FluidSimulation_ObjectInteraction(Scene* scene) {
	real h, dt, cfl, cfl_interface;
	real Lx, Ly, Lz;
	int nX, nY, nZ;
	int endFrame;
	real narrowband_srl = 4;

	const unsigned int redist_freq = 7;
	const unsigned int redist_min = 4, redist_max = 25;

	endFrame = 100;
	dt = 1.0f/30.0f;	
	cfl = 3.5f; cfl_interface = 2.0f;

	nX = nZ = 16 / 1; nY = 16 / 1;
	h = 1.0 / double(nX);	

	Lx = nX*h; Ly = nY*h; Lz = nZ*h;

	FluidScene *fluidscene = new FluidScene(dt);
	fluidscene->initFluidSceneIntegrator(new FluidSceneIntegrator_StaticRigidObject(dt));	
	fluidscene->setEndFrame(endFrame);
	//fluidscene->enableSaveBinary(500);
	fluidscene->enableSaveRenderObjectsOnly(500);
	//fluidscene->enableSaveImage();		
	scene->push_SceneObject(fluidscene);

	Liquid3D *liquid3D = new Liquid3D(nX, nY, nZ, h, dt, cfl, cfl_interface);
	liquid3D->release();

	////////////////////////////////////////////////////////////////////////	
	// Interface Tracking
	liquid3D->tracking = new InterfaceTracking_SRL3D(nX, nY, nZ, h, narrowband_srl);

	////////////////////////////////////////////////////////////////////////
	// CFD Integrator
	CFDIntegrator_Liquid* itg = new CFDIntegrator_Liquid();	

	itg->set_CFDAddForce(new CFDAddForce_Gravity(-9.81));
	itg->set_CFDAdvection(new CFDAdvection_SemiLagrangian<>());
	
	CFDProjection_Liquid* project = new CFDProjection_FreeSurface_Gibou();
	
	project->set_InterfaceTracking(liquid3D->tracking);
	project->init_PoissonSolver_GMG(liquid3D->solver, 12); project->set_iteration_number_GMG(3,5);
	project->set_divergence_threshold(1.0e-3);

	itg->set_CFDProjection(project);
	itg->set_CFDTimeStepRestriction(new CFDTimeStepRestriction(cfl));	

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

	//LiquidSource_Continuous* source_c = new LiquidSource_Continuous(new SphereDistanceField3D<real>(Vector3<real>(0.5f*Lx, 0.83f*Ly, 0.5f*Lz), 0.05f*Lx));
	//source_c->add_interval(1,3);
	//source_c->set_vel(Vector2<real>(0,-2));
	//source = source_c;

	//LiquidSource_Discrete* source_d = new LiquidSource_Discrete(new SphereDistanceField3D<real>(Vector3<real>(0.5f*Lx, 0.83f*Ly, 0.5f*Lz), 0.10f*Lx));
	//source_d->add_time_point(0.5); source_d->add_time_point(1.0); source_d->add_time_point(1.5);
	//source = source_d;

	itg->set(liquid3D->tracking, itg_srl, extrapolation, source);	
	liquid3D->integrator = itg;

	////////////////////////////////////////////////////////////////////////	
	// Set Phi
	liquid3D->setPhi(WaterDropDistanceField3D<real>(Vector3<real>(0.5f*Lx, 0.83*Ly, 0.5f*Lz), 0.15f*Ly, 0.10f*Ly));

	//GridDistanceField3D<real> grid_field; 
	//grid_field.set("../Data/bunny.obj", 200, 200, 200, "../Data/bunny_phi");
	//grid_field.set("../Data/Armadillo.ply", 200, 200, 200, "../Data/Armadillo_phi", true, arithmetic::rotateY(real(180.0)));
	//grid_field.load("../Data/Armadillo_phi"); grid_field.set_rigid_transformation(translate(Vector3<real>(0.5,0.5,0.5)));

	//grid_field.load("../Data/bunny_phi"); 
	//grid_field.set_transformation(geometryspace::translate(Vector3<real>(0.5,0.5,0.5)));
	//liquid3D->setPhi(grid_field);

	fluidscene->initFluidSimulator(liquid3D);	

	////////////////////////////////////////////////////////////////////////	
	// Set Object
	//FluidScene_Objects3D* object3D = new FluidScene_Objects3D();
	//object3D->push_DynamicObject(StaticRigidBody3D<T>(Sphere<T>(Vector3<T>(nX*h*0.5f, nY*h*0.42f, nZ*h*0.5f), nY*h*0.10f)));
	//fluidscene->initFluidScene_Objects(object3D);
}

void SRL_FluidSimulation_SmallScale(Scene* scene) {

	const real g = - 9.81, sigma = 0.0005;

	const int endFrame = 400;
	const real dt = 1.0 / 300.0;	
	const real cfl = 1.0f, cfl_interface = 1.0f;

	const real narrowband_srl = 4;
	const unsigned int redist_freq = 8;
	const unsigned int redist_min = 4, redist_max = 50;

	const int ratio = 1;
	const int nX = 64*ratio, nY = 64*ratio, nZ = nX; 
	const real h = 0.03 / double(nX);	

	const real Lx = nX*h, Ly = nY*h, Lz = nZ*h;

	FluidScene *fluidscene = new FluidScene(dt);
	fluidscene->initFluidSceneIntegrator(new FluidSceneIntegrator_StaticRigidObject(dt));	
	fluidscene->setEndFrame(endFrame);
	fluidscene->enableSaveBinary(500);
	fluidscene->enableSaveRenderObjectsOnly(500);
	fluidscene->enableSaveImage();		
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

	CFDProjection_Liquid* project = new CFDProjection_FreeSurface_Gibou(true, sigma);
	project->set_InterfaceTracking(liquid3D->tracking);
	project->init_PoissonSolver_GMG(liquid3D->solver, 12); project->set_iteration_number_GMG(3,5);
	project->set_divergence_threshold(1.0e-6);
	itg->set_CFDProjection(project);

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
	
	Levelset_Extrapolation<real>* extrapolation = NULL;
	extrapolation = new Levelset_Extrapolation_PDE<real>(narrowband_srl, 0.5f);

	itg->set(liquid3D->tracking, itg_srl, extrapolation);	
	liquid3D->integrator = itg;

	////////////////////////////////////////////////////////////////////////	
	// Set Phi
	liquid3D->setPhi(WaterDropDistanceField3D<real>(Vector3<real>(0.5f*Lx, 0.83*Ly, 0.5f*Lz), 0.13f*Lx, 0.05f*Ly));
	fluidscene->initFluidSimulator(liquid3D);	
}
