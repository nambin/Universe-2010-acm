#include <stdafx.h>
#include "SRL_Experiment.h"

using namespace std;
using namespace geometryspace;
using namespace rigidspace;
using namespace fluidspace;

typedef fluidspace::real T;

void SRL_FluidSimulation_3D_Demo_03(Scene* scene) {
#ifdef USE_OPENMP	
	const int maxthreads = omp_get_max_threads();
	const int numThreads = maxthreads;

	omp_set_num_threads(numThreads);
	log_info << endl << "OpenMP " << numThreads << " Threads are initiated" << endl << endl;
#endif

	const real g = - 9.81, sigma = 0.000;

	const int endFrame = 300;
	const real dt = 1.0 / 30.0;	
	const real cfl = 3.5, cfl_interface = 1.5;

	const real narrowband_srl = 4;
	const unsigned int redist_freq = 15;
	const unsigned int redist_min = 5, redist_max = 50;

	const real ratio = 2*2*2*2*2;
//	const real ratio = 2*2*2;
	const int nX = 7 * ratio, nY = 7 * ratio, nZ = 5 * ratio;
	const real h = 4.5 / double(nX);	

	const real Lx = nX*h, Ly = nY*h, Lz = nZ*h;

	FluidScene *fluidscene = new FluidScene(dt);
	fluidscene->initFluidSceneIntegrator(new FluidSceneIntegrator_StaticRigidObject(dt));	
	fluidscene->setEndFrame(endFrame);
	fluidscene->enableSaveBinary(500);
	fluidscene->enableSaveRenderObjectsOnly(500);
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

	itg->set_CFDAddForce(new CFDAddForce_Gravity(-9.81));
	itg->set_CFDAdvection(new CFDAdvection_SemiLagrangian<>());

	CFDProjection_Liquid* project = new CFDProjection_FreeSurface_Gibou(false, sigma);

	project->set_InterfaceTracking(liquid3D->tracking);
	project->init_PoissonSolver_GMG(liquid3D->solver, 12); project->set_iteration_number_GMG(3,5);
	project->set_divergence_threshold(1.0e-3);

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

	Levelset_Extrapolation<real>* extrapolation = new Levelset_Extrapolation_PDE<real>(narrowband_srl, 0.5f);

	LiquidSource* source = NULL;

	itg->set(liquid3D->tracking, itg_srl, extrapolation, source);	
	liquid3D->integrator = itg;

	const real r_dam_pos_x = 0.2f, r_dam_width = 0.1f;
	const real dam_pos_x = (r_dam_pos_x * nX) * h;

	////////////////////////////////////////////////////////////////////////	
	// Set Phi	
	liquid3D->setPhi(CubeDistanceField3D<real>(Vector3<real>(-10.0f * Lx), Vector3<real>(dam_pos_x, 0.9 * (real(48)/real(64)) * Ly, 10.0f*Lz)));
	fluidscene->initFluidSimulator(liquid3D);	

	////////////////////////////////////////////////////////////////////////	
	// Set Object
	FluidScene_Objects3D* object3D = new FluidScene_Objects3D();

	const T r_dam_base_height = 0.1f * (real(58)/real(64));
	const T r_dam_pos_y = 0.45f * (real(50)/real(64)), r_dam_height = 0.95 - r_dam_pos_y;
	const T r_dam_base_z = 0.0f;

	const T dam_width = int(r_dam_width * nX) * h, dam_pos_y = int(r_dam_pos_y * nY) * h;
	const T dam_base_height = int(r_dam_base_height * nY) * h, dam_height = int(r_dam_height * nY) * h;
	const T dam_base_z = int(r_dam_base_z * nZ) * h;

	Vector3<T> min_p = Vector3<T>(dam_pos_x,0,dam_base_z), max_p = Vector3<T>(dam_pos_x + dam_width, dam_base_height, Lz-dam_base_z);
	Cube<T> dam_base(min_p, max_p);
	object3D->push_DynamicObject(StaticRigidBody3D<T>(dam_base));

	log_info << "Dam Base _ " << min_p << " " << max_p << std::endl;
	log_info << "Dam Base _ " << (min_p+max_p) * real(0.5) << " " << max_p - min_p << std::endl;

	min_p = Vector3<T>(dam_pos_x,dam_pos_y,0), max_p = Vector3<T>(dam_pos_x + dam_width, Ly, Lz);
	Cube<T> dam(min_p, max_p);
	object3D->push_DynamicObject(StaticRigidBody3D<T>(dam));

	log_info << "Dam Upper _ " << min_p << " " << max_p << std::endl;
	log_info << "Dam Upper _ " << (min_p+max_p) * real(0.5) << " " << max_p - min_p << std::endl << std::endl;

	unsigned int num_pole_x = 3;
	unsigned int num_pole_z = 4;

	const real pole_radius = Lz * 0.05;
	const real pole_length_min = Ly * 0.2 * (real(62)/real(64));
	const real pole_length_max = Ly * 0.3 * (real(62)/real(64));

	const real _xs = dam_pos_x + dam_width;
	const real _xe = Lx;

	for(unsigned x=0;x<num_pole_x;++x) {
		for(unsigned z=0;z<num_pole_z;++z) {
			Vector3<real> pole_center;

			real xs, xe, zs, ze;

			if(x == 1 && z % 2 == 0)
				continue;

			if( x == 0 ) {
				xs = _xs + (_xe-_xs) * real(x+0) / real(num_pole_x) + pole_radius * 3.0;
				xe = _xs + (_xe-_xs) * real(x+1) / real(num_pole_x) - pole_radius * 1.0;

				pole_center[0] = xs + (xe-xs) * random();
			}
			else if(x == 1) {
				xs = _xs + (_xe-_xs) * real(x+0) / real(num_pole_x) + pole_radius * 1.0;
				xe = _xs + (_xe-_xs) * real(x+1) / real(num_pole_x) - pole_radius * 1.0;

				pole_center[0] = xs + (xe-xs) * 0.5;
			}
			else if(x == 2) {
				xs = _xs + (_xe-_xs) * real(x+0) / real(num_pole_x) + pole_radius * 2.0;
				xe = _xs + (_xe-_xs) * real(x+1) / real(num_pole_x) - pole_radius * 2.0;

				pole_center[0] = xs + (xe-xs) * random();
			}			
			
			if(x == 0 || x == 2) {
				zs = Lz * real(z+0) / real(num_pole_z) + pole_radius * 1.5;
				ze = Lz * real(z+1) / real(num_pole_z) - pole_radius * 1.5;

				pole_center[2] = zs + (ze-zs) * random();
			}
			else {
				zs = Lz * real(z/2+0) / real(2);
				ze = Lz * real(z/2+1) / real(2);

				pole_center[2] = zs + (ze-zs) * 0.5;
			}						

			real pole_length = pole_length_min + (pole_length_max - pole_length_min) * random();

			if(x == 0)		pole_length /= 1.5;
			else if(x == 1)	pole_length /= 1.3;

			pole_center[1] = pole_length * 0.5;

			Cylinder<T> pole(pole_radius, pole_length);
			pole.set_rigid_transformation(arithmetic::translate(pole_center));		

			log_info << "Cylinder " << " _ " << "r _ " << pole_radius << " , length _ " << pole_length << std::endl;
			log_info << "Translate _ " << pole_center << " Rotation _ 0" << std::endl << std::endl;

			object3D->push_DynamicObject(StaticRigidBody3D<T>(pole));
		}
	}
	
	fluidscene->initFluidScene_Objects(object3D);
}