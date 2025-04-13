#include <stdafx.h>
#include "SRL_Experiment.h"

using namespace std;
using namespace geometryspace;
using namespace rigidspace;
using namespace fluidspace;

typedef fluidspace::real T;

class LiquidSource_Falling_Star : public LiquidSource {
public : 
	LiquidSource_Falling_Star() {
		star.load("../Data/polyhedra_star_phi"); 
	}		

	void set_vertical_falling(real Lx, real Ly, real Lz, real ts, real te, real scale_min, real scale_max, real v_min, real v_max, int num_per_sec) {
		const int num = (te-ts) * num_per_sec;

		t0_position.resize(num);
		t0_scale.resize(num);
		t0_time_point.resize(num);
		vel0.resize(num); 

		for(unsigned int i=0;i<num;++i) {
			const real t_interval = (te - ts) / real(num);

			if(i==0) t0_time_point[i] = ts + 0.001;
			else	 t0_time_point[i] = (random.getFloat() * t_interval) + (i * t_interval) + ts;

			t0_scale[i] = (random.getFloat() * (scale_max - scale_min) + scale_min) * MIN(Lx, Lz);
			
			vel0[i][1] = (random.getFloat() * (v_max - v_min) + v_min); 
			vel0[i][0] = vel0[i][2] = 0;

			bool not_overlap = false;

			const real eps = t0_scale[i] * 0.6;			

			const real x_max = Lx - eps, x_min = eps;
			const real z_max = Lz - eps, z_min = eps;

			while(not_overlap == false) {

				t0_position[i][0] = (random.getFloat() * (x_max - x_min) + x_min);
				t0_position[i][2] = (random.getFloat() * (z_max - z_min) + z_min);
				t0_position[i][1] = Ly - eps;

				if( i != 0 && ::abs(t0_position[i][0] - t0_position[i-1][0]) < eps  && ::abs(t0_position[i][2] - t0_position[i-1][2]) < eps)	not_overlap = false;
				else																															not_overlap = true;

			}
		}
	}

	void add(InterfaceTracking* tracking, CFDsolver& cfd, real dt) const {

		for(size_t i=0;i<t0_time_point.size();i++) {
			if(time <= t0_time_point[i] && t0_time_point[i] < time + dt) {								
				star.set_rigid_transformation(arithmetic::translate(t0_position[i]));
				star.set_scale(t0_scale[i]);
				tracking->addSource(star);

				CFDsolver3D& cfd3D = (CFDsolver3D&)cfd;
				setVel(cfd3D, vel0[i], star);		
			}
		}

		time += dt;
	}	

protected : 
	basis::RandomGen<real>				random;

	std::vector<Vector3<real>>			t0_position;
	std::vector<real>					t0_scale;
	std::vector<real>					t0_time_point;
	std::vector<Vector3<real>>			vel0;

	mutable GridDistanceField3D<real>	star;
};

void SRL_FluidSimulation_3D_Demo_04(Scene* scene) {
#ifdef USE_OPENMP	
	const int maxthreads = omp_get_max_threads();
	const int numThreads = maxthreads;

	omp_set_num_threads(numThreads);
	log_info << endl << "OpenMP " << numThreads << " Threads are initiated" << endl << endl;
#endif

	const real g = - 9.81, sigma = 0.000;

	const int endFrame = 360;
	const real dt = 1.0 / 60.0;	
	const real cfl = 3.5, cfl_interface = 1.5;

	const real narrowband_srl = 4;
	const unsigned int redist_freq = 15;
	const unsigned int redist_min = 5, redist_max = 50;

	const real ratio = 2*2*2*2*2;
//	const real ratio = 2*2;
	const int nX = 6 * ratio, nY = 6 * ratio, nZ = 4 * ratio;
	const real h = 2.5 / double(nX);	

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

	LiquidSource_Falling_Star* source = new LiquidSource_Falling_Star();
	source->set_vertical_falling(Lx, Ly, Lz, 0, 9, 0.3, 0.4, -4.0, -5.0, 8);
	//source->set_vertical_falling(Lx, Ly, 0.01, 5.5, 0.5, 0.5, -4.0, 5);

	itg->set(liquid3D->tracking, itg_srl, extrapolation, source);	
	liquid3D->integrator = itg;

	////////////////////////////////////////////////////////////////////////	
	// Set Phi
	liquid3D->setPhi(WaterDropDistanceField3D<real>(Vector3<real>(10000.0), 0.0, 0.04*Ly));

	fluidscene->initFluidSimulator(liquid3D);	
}