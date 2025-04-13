#include <stdafx.h>
#include "SRL_Experiment.h"

using namespace std;
using namespace geometryspace;
using namespace rigidspace;
using namespace fluidspace;

typedef fluidspace::real T;

void SRL_FluidSimulation_ObjectInteraction2D(Scene* scene) {

	const real g = - 9.81;

	const int endFrame = 200;
	const real dt = 1.0 / 30.0;	
	const real cfl = 3.1, cfl_interface = 1.6;

	const real narrowband_srl = 4;
	const unsigned int redist_freq = 7;
	const unsigned int redist_min = 4, redist_max = 25;

	const int ratio = 1;
	const int nX = 64 * ratio, nY = 96 * ratio;
	const real h = 1.0 / double(nX);	

	const real Lx = nX*h, Ly = nY*h;

	FluidScene *fluidscene = new FluidScene(dt);
	fluidscene->initFluidSceneIntegrator(new FluidSceneIntegrator_StaticRigidObject(dt));	
	fluidscene->setEndFrame(endFrame);
	fluidscene->enableSaveImage();		
	scene->push_SceneObject(fluidscene);

	Liquid2D *liquid2D = new Liquid2D(nX, nY, h, dt, cfl, cfl_interface);
	liquid2D->release();

	////////////////////////////////////////////////////////////////////////	
	// Interface Tracking
	liquid2D->tracking = new InterfaceTracking_SRL2D(nX, nY, h, narrowband_srl);

	////////////////////////////////////////////////////////////////////////
	// CFD Integrator
	CFDIntegrator_Liquid* itg = new CFDIntegrator_Liquid();	

	itg->set_CFDAddForce(new CFDAddForce_Gravity(-9.81));
	itg->set_CFDAdvection(new CFDAdvection_SemiLagrangian<>());

	CFDProjection_Liquid* project = new CFDProjection_FreeSurface_Gibou();
//	project->set_liquid_biased_filter(true);

	project->set_InterfaceTracking(liquid2D->tracking);
	project->init_PoissonSolver_GMG(liquid2D->solver, 12); project->set_iteration_number_GMG(3,5);

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

	LiquidSource_Continuous* source_c = new LiquidSource_Continuous(new CircleDistanceField2D<real>(Vector2<real>(0.5f*Lx, 0.8f*Ly), 0.05f*Lx));
	source_c->add_interval(1,3);
	source_c->set_vel(Vector2<real>(0,-2));
	source = source_c;

	//LiquidSource_Discrete* source_d = new LiquidSource_Discrete(new CircleDistanceField2D<real>(Vector2<real>(0.5f*Lx, 0.8f*Ly), 0.15f*Lx));
	//source_d->add_time_point(0.5); source_d->add_time_point(1.0); source_d->add_time_point(1.5);
	//source = source_d;

	itg->set(liquid2D->tracking, itg_srl, extrapolation, source);	
	liquid2D->integrator = itg;

	////////////////////////////////////////////////////////////////////////	
	// Set Phi
	liquid2D->setPhi(WaterDropDistanceField2D<real>(Vector2<real>(0.5f*Lx, 0.80*Ly), 0.18f*Lx, 0.10f*Ly));
	fluidscene->initFluidSimulator(liquid2D);	

	////////////////////////////////////////////////////////////////////////	
	// Set Object
	FluidScene_Objects2D* object2D = new FluidScene_Objects2D();
	//object2D->push_DynamicObject(StaticRigidBody2D<T>(Circle<T>(Vector2<T>(Lx*0.5f, Ly*0.42f), Lx*0.12f)));
	//object2D->push_DynamicObject(StaticRigidBody2D<T>(Circle<T>(Vector2<T>(Lx*0.1f, Ly*0.1f), Lx*0.15f)));
	object2D->push_DynamicObject(StaticRigidBody2D<T>(Circle<T>(Vector2<T>(Lx*0.9f, Ly*0.1f), Lx*0.15f)));
	fluidscene->initFluidScene_Objects(object2D);
}

void SRL_FluidSimulation_SmallScale2D(Scene* scene) {

	real g = -9.81, sigma = 0.0728;
	real rho_neg = 1000, rho_pos = 1.226;

	const int endFrame = 250;
	const real dt = 1.0f/500.0f;	
	const real cfl = 10.0f, cfl_interface = 1.0f;

	const real narrowband = 4;
	const unsigned int redist_freq = 12;
	const unsigned int redist_min = 4, redist_max = 50;

	const int ratio = 1;
	const int nX = 64*ratio, nY = 96*ratio; 
	const real h = 0.05 / double(nX);	

	const real Lx = nX*h, Ly = nY*h;

	FluidScene *fluidscene = new FluidScene(dt);
	fluidscene->initFluidSceneIntegrator(new FluidSceneIntegrator_StaticRigidObject(dt));	
	fluidscene->setEndFrame(endFrame);
	//	fluidscene->enableSaveBinary(500);
	//	fluidscene->enableSaveRenderObjectsOnly(500);
	fluidscene->enableSaveImage();		
	scene->push_SceneObject(fluidscene);

	Liquid2D *liquid2D = new Liquid2D(nX, nY, h, dt, cfl, cfl_interface);
	liquid2D->release();

	////////////////////////////////////////////////////////////////////////	
	// Interface Tracking
	liquid2D->tracking = new InterfaceTracking_SRL2D(nX, nY, h, narrowband);

	////////////////////////////////////////////////////////////////////////
	// CFD Integrator
	CFDIntegrator_Liquid* itg = new CFDIntegrator_Liquid();	

	//itg->set_CFDAddForce(new CFDAddForce_Gravity(g));
	itg->set_CFDAddForce(new CFDAddForce_SurfaceTension(g, sigma, rho_neg, rho_pos));
	itg->set_CFDAdvection(new CFDAdvection_SemiLagrangian<>());

	//CFDProjection_Liquid* project = new CFDProjection_FreeSurface_Gibou(false, sigma);
	CFDProjection_Liquid* project = new CFDProjection_CSF(rho_neg, rho_pos);	
	project->set_InterfaceTracking(liquid2D->tracking);
	project->init_PoissonSolver_GMG(liquid2D->solver, 12); project->set_iteration_number_GMG(3,5);
	project->set_divergence_threshold(1.0e-6);
	itg->set_CFDProjection(project);

	itg->set_CFDTimeStepRestriction(new CFDTimeStepRestriction_SurfaceTension(cfl, g, sigma));
	itg->set_CFDVelocityBoundaryCondition(new GridBoundaryCondition_U_FreeSlip(), new GridBoundaryCondition_V_FreeSlip(), new GridBoundaryCondition_W_FreeSlip());

	////////////////////////////////////////////////////////////////////////	
	// Interface Integrator
	InterfaceIntegrator_SRL* itg_srl = new InterfaceIntegrator_SRL();

	SRL_Redistancing* redist = new SRL_Redistancing_PDE_WENO(narrowband_srl, 0.5, true, 0);
	redist->set_redistance_cfl_frequency(redist_freq);
	redist->set_redistance_min_max_updates(redist_min,redist_max);

	itg_srl->set_Advection(new CFDAdvection_SemiLagrangian_SRL<CFDTracking_RK2>(true));		
	itg_srl->set_Redistancing(redist);		
	itg_srl->set_TimeStepRestriction(new CFDTimeStepRestriction(cfl_interface));

	Levelset_Extrapolation<real>* extrapolation = NULL;

	if(dynamic_cast<CFDProjection_CSF*>(project) == NULL)
		extrapolation = new Levelset_Extrapolation_PDE<real>(narrowband, 0.5f);

	itg->set(liquid2D->tracking, itg_srl, extrapolation);	
	liquid2D->integrator = itg;

	////////////////////////////////////////////////////////////////////////	
	// Set Phi
	liquid2D->tracking->setPhi(WaterDropDistanceField2D<real>(Vector2<real>(0.5f*Lx, 0.6f*Ly), 0.12f*Ly, 0.2f*Ly));
	liquid2D->tracking->createDrawableInterface();
	
	fluidscene->initFluidSimulator(liquid2D);	
}
