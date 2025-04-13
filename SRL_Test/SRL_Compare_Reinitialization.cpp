#include <stdafx.h>
#include "SRL_Experiment.h"

using namespace std;
using namespace geometryspace;
using namespace rigidspace;
using namespace fluidspace;

typedef fluidspace::real T;

void SRL_Compare_Reinitialization_RotatingBunny(Scene* scene) {

#ifdef USE_OPENMP	
	int maxthreads = omp_get_max_threads();
	int numThreads = maxthreads;

	omp_set_num_threads(numThreads);
	log_info << endl << "OpenMP " << numThreads << " Threads are initiated" << endl << endl;
#endif

	const unsigned int narrowband_srl = 5;

	const unsigned int n_rotation = 1;

	const unsigned int n_step = 950;
	const unsigned int res = 128;
	const unsigned int endTime = 200; 
	const T dt = float(endTime) / (T(n_step) * 1.0);

	const unsigned int endframe = int(endTime / dt + 1) * n_rotation + n_rotation + 1;
	
	const grid_index_type nX(res), nY(res), nZ(res);
	const real L = 100;
	const real h = L / real(res);

	InterfaceTrackingTest test;
	std::string scheme;

	for(int i=0;i<3;++i) {

		test.release();

		if(i==0)	scheme = "WENO3rd";
		if(i==1)	scheme = "Direct";
		if(i==2)	scheme = "WENO";		

		log_info.set_filename("InterfaceTracking_Rotating Bunny_" + scheme);	test.test_name = "RotatingBunny" + scheme;

		////////////////////////////////////////////////////////////////////////
		// Geometry & Velocity Field
		ZalesakDistanceField3D<real> zalesak(L);	
		test.vel = new ZalesakVelocityField3D<real>(endTime, L);

		//////////////////////////////////////////////////////////////////////////
		// Frame
		test.FluidSimulator::set(dt, 1.0e+6);
		test.FluidSimulator::setEndFrame(endframe);	
		test.saveBinaries = true;
		test.saveImages = false;

		GridDistanceField3D<real> bunny; bunny.load("../Data/bunny_phi"); 
		bunny.set_scale(L*0.45);
		bunny.set_rigid_transformation(arithmetic::translate(Vector3<real>(0.5*L,0.75*L,0.5*L)));

		InterfaceTracking_SRL3D* srl3D = new InterfaceTracking_SRL3D();
		srl3D->set(nX, nY, nZ, h, narrowband_srl, false);
		srl3D->setPhi(bunny);
		test.tracking = srl3D;
		
		InterfaceIntegrator_SRL* itgSRL = new InterfaceIntegrator_SRL();
		itgSRL->set_Advection(new CFDAdvection_SemiLagrangian_SRL<CFDTracking_RK4>(true));

		if( i == 0 )		itgSRL->set_Redistancing(new SRL_Redistancing_PDE_WENO3rd(narrowband_srl, 0.4, true, 0));
		else if( i == 1 )	itgSRL->set_Redistancing(new SRL_Redistancing_Direct(narrowband_srl));
		else if( i == 2 )	itgSRL->set_Redistancing(new SRL_Redistancing_PDE_WENO(narrowband_srl, 0.4, true, 0));		

		itgSRL->set_TimeStepRestriction(new CFDTimeStepRestriction(1.0e+6));
		test.integrator = itgSRL;

		while(test.frame < test.endFrame) {
			if(test.frame % 100 == 0 || (test.frame % (n_step / 4)) == 0 || test.frame == 0 || (test.frame % n_step == 0)) {
				test.tracking->createDrawableInterface(Grid3b());

				char str[128];
				sprintf_s(str, 128, "%04d", test.frame);
				test.tracking->exportOBJ(output_folder + "OBJs/" + "Bunny" + "_" + scheme + "_" + str + ".obj");
			}				

			test.updateOneStep();
		}

		//while(test.frame < test.endFrame) 
		//	test.updateOneStep();
	}	
}

void SRL_Compare_Reinitialization_RotatingZalesak(Scene* scene) {

#ifdef USE_OPENMP	
	int maxthreads = omp_get_max_threads();
	int numThreads = maxthreads;

	omp_set_num_threads(numThreads);
	log_info << endl << "OpenMP " << numThreads << " Threads are initiated" << endl << endl;
#endif

	const unsigned int narrowband_srl = 5;

	const unsigned int n_rotation = 1;

	const unsigned int res = 100;
	const unsigned int endTime = 200; 
	const unsigned int n_step = 850;
	const T dt = 200.0 / n_step;

	const unsigned int endframe = int(endTime / dt + 1) * n_rotation + n_rotation + 1;
	
	const grid_index_type nX(res), nY(res), nZ(res);
	const real L = 100;
	const real h = L / real(res);

	InterfaceTrackingTest test;

	for(int i=0;i<2;++i) {

		if(i == 1 || i == 2) continue;

		test.release();

		string scheme;
		switch(i) {
			case 0 : scheme = "WENO3rd";	break;	
			case 1 : scheme = "WENO";	break;	
			case 2 : scheme = "Direct";	break; 
		}

		if( i == 0)		 { log_info.set_filename("InterfaceTracking _ Zalesak _ WENO 3rd");	test.test_name = "ZalesakWENO3rd";	}
		else if( i == 1) { log_info.set_filename("InterfaceTracking _ Zalesak _ WENO");	test.test_name = "ZalesakWENO";			}
		else			 { log_info.set_filename("InterfaceTracking _ Zalesak _ Direct"); test.test_name = "ZalesakDirect";		}

		////////////////////////////////////////////////////////////////////////
		// Geometry & Velocity Field
		ZalesakDistanceField3D<real> zalesak(L);	
		test.vel = new ZalesakVelocityField3D<real>(endTime, L);

		//////////////////////////////////////////////////////////////////////////
		// Frame
		test.FluidSimulator::set(dt, 1.0e+6);
		test.FluidSimulator::setEndFrame(endframe);	
		test.saveBinaries = true;
		test.saveImages = false;

		InterfaceTracking_SRL3D* srl3D = new InterfaceTracking_SRL3D();
		srl3D->set(nX, nY, nZ, h, narrowband_srl, false);
		srl3D->setPhi(zalesak);
		test.tracking = srl3D;

		InterfaceIntegrator_SRL* itgSRL = new InterfaceIntegrator_SRL();
		itgSRL->set_Advection(new CFDAdvection_SemiLagrangian_SRL<CFDTracking_RK4>(true));

		if( i == 0)			itgSRL->set_Redistancing(new SRL_Redistancing_PDE_WENO3rd(narrowband_srl, 0.3, true, 0));
		else if( i == 1)	itgSRL->set_Redistancing(new SRL_Redistancing_PDE_WENO(narrowband_srl, 0.3, true, 0));
		else				itgSRL->set_Redistancing(new SRL_Redistancing_Direct(narrowband_srl));

		itgSRL->set_TimeStepRestriction(new CFDTimeStepRestriction(1.0e+6));
		test.integrator = itgSRL;

		while(test.frame < test.endFrame) {
			if(test.frame % 50 == 0 || (test.frame % (n_step / 4)) == 0) {
				test.tracking->createDrawableInterface(Grid3b());

				char str[128];
				sprintf_s(str, 128, "%04d", test.frame);
				test.tracking->exportOBJ(output_folder + "OBJs/" + "Zalesak" + "_" + scheme + "_" + str + ".obj");
			}				

			test.updateOneStep();
		}			
	}
}

void SRL_Compare_Reinitialization_2DTest(Scene* scene) {

#ifdef USE_OPENMP	
	int maxthreads = omp_get_max_threads();
	int numThreads = maxthreads;

	omp_set_num_threads(numThreads);
	log_info << endl << "OpenMP " << numThreads << " Threads are initiated" << endl << endl;
#endif

	const unsigned int narrowband_srl = 5;
	const unsigned int n_rotation = 5;

	//const unsigned int res = 80;	const unsigned int n_step = 628;
	const unsigned int res = 101; const unsigned int n_step = 628;

	const unsigned int endTime = 1.0f; 
	const T dt = T(endTime) / T(n_step);

	const unsigned int endframe = int(endTime / dt + 1) * n_rotation + n_rotation + 1;

	const grid_index_type nX(res), nY(res), nZ(res);
	const real L = 100;
	const real h = L / real(res);

	InterfaceTrackingTest* test = new InterfaceTrackingTest();

	for(int i=0;i<4;++i) {	
//		DistanceField2D<real> * distfield = new CircleDistanceField2D<real>(Vector2<real>(0.5f*L, 0.75f*L), 0.12f*L);
		DistanceField2D<real> * distfield = new ZalesakDistanceField2D<real>(L);
		VelocityField2D<real> * velfield = new RotationVelocityField2D<real>(endTime, L, L);

		test->release();

		string scheme;
		switch(i) {
			case 0 : scheme = "Direct";		break; 
			case 1 : scheme = "WENO";		break;				
			case 2 : scheme = "WENO3rd";	break;	
			case 3 : scheme = "NoRedist";	break;
		}

		log_info.set_filename("InterfaceTracking_Zalesak2D_" + scheme);	
		test->test_name = "Zalesak2D_" + scheme;		

		////////////////////////////////////////////////////////////////////////
		// Geometry & Velocity Field		
		test->vel = velfield;

		//////////////////////////////////////////////////////////////////////////
		// Frame
		test->FluidSimulator::set(dt, 1.0e+6);
		test->FluidSimulator::setEndFrame(endframe);	
		test->saveBinaries = false;
		test->saveImages = false;

		InterfaceTracking_SRL2D* srl2D = new InterfaceTracking_SRL2D();
		srl2D->set(nX, nY, h, narrowband_srl, false);
		srl2D->setPhi(*distfield);
		test->tracking = srl2D;

		InterfaceIntegrator_SRL* itgSRL = new InterfaceIntegrator_SRL();
		itgSRL->set_Advection(new CFDAdvection_SemiLagrangian_SRL<CFDTracking_RK4>(true));

		if( i == 0 )		itgSRL->set_Redistancing(new SRL_Redistancing_Direct(narrowband_srl));
		else if( i == 1 )	itgSRL->set_Redistancing(new SRL_Redistancing_PDE_WENO(narrowband_srl, 0.4, true, 0));
		else if( i == 2 )	itgSRL->set_Redistancing(new SRL_Redistancing_PDE_WENO3rd(narrowband_srl, 0.4, true, 0));
		else				itgSRL->set_Redistancing(new SRL_Redistancing_No());

		itgSRL->set_TimeStepRestriction(new CFDTimeStepRestriction(1.0e+6));
		test->integrator = itgSRL;
		test->tracking->createDrawableInterface(Grid2b());		

		while(test->frame < test->endFrame) {
			if(test->frame == 0 || (test->frame % n_step) == 0) {
				char frame_s[64];
				sprintf_s(frame_s, "%03d", test->frame);

				char res_s[64];
				sprintf_s(res_s, "%03d", res);

				if(test->frame == 0)				srl2D->srl->phi.save_as_matlab((output_folder + "Matlab/redist_" + res_s + "_Exact.m").c_str(), "Exact");
				else if(test->frame % n_step == 0) 	srl2D->srl->phi.save_as_matlab((output_folder + "Matlab/redist_" + res_s + "_" + frame_s + "_" + scheme + ".m").c_str(), scheme.c_str());
			}	

			test->updateOneStep();					
		}
	}

	scene->push_SceneObject(test);
	test->enableSaveImage();
}