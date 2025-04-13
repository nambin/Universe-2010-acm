#include <stdafx.h>
#include "SRL_Experiment.h"

using namespace std;
using namespace geometryspace;
using namespace rigidspace;
using namespace fluidspace;

typedef fluidspace::real T;

void SRL_Compare_AdvectionTest_Rotation(Scene* scene) {

#ifdef USE_OPENMP	
	int maxthreads = omp_get_max_threads();
	int numThreads = maxthreads;

	omp_set_num_threads(numThreads);
	log_info << endl << "OpenMP " << numThreads << " Threads are initiated" << endl << endl;
#endif

	const unsigned int res_basic = 100;
	const unsigned int n_step_basic = 628;

	const unsigned int narrowband_srl = 4;
	const unsigned int n_rotation = 1;
	
	const T endTime = 1.0, L = 100.0; 

	InterfaceTrackingTest* test = NULL;

	for(int i=0;i<5;++i) {		
		for(int j=0;j<2;++j) {

			if( i == 0 || i == 1 || i == 2 || i == 4 || j == 1 )
				continue;

			test = new InterfaceTrackingTest();

			string scheme;
			string distancefield;

			switch(i) {
				case 0 : scheme = "SL";		break;	case 1 : scheme = "MacCormack";	break; 
				case 2 : scheme = "BFECC";	break;
				case 3 : scheme = "SLC";	break;	case 4 : scheme = "SRL";	break;
			}
			
			unsigned int res, n_step;
			switch(i) {
				case 0 : 
				case 1 : 
				case 2 : 
				case 3 : res = res_basic * 4;	n_step = n_step_basic;	break;
				case 4 : res = res_basic;		n_step = n_step_basic;	break;
			}

			bool need_limiting;
			switch(j) {
				case 0 : distancefield = "Zalesak";		need_limiting = false;	break;
				case 1 : distancefield = "Bunny";		need_limiting = true;	break;
			}

			log_info.set_filename(distancefield + "_" + scheme);
			test->test_name = distancefield + "_" + scheme;
			test->image_filename = distancefield + "_" + scheme;

			const T dt = real(endTime) / real(n_step);
			const unsigned int endframe = n_step_basic * n_rotation + 2;
			const real h = L / real(res);

			log_info << "Number of Steps _ " << n_step << endl;
			log_info << "Number of Rots  _ " << n_rotation << endl;
			log_info << "Time Step       _ " << dt << endl;

			//////////////////////////////////////////////////////////////////////////
			// Frame
			test->FluidSimulator::set(dt, 1.0e+6);
			test->FluidSimulator::setEndFrame(endframe);	

			test->saveBinaries = test->saveImages = false;

			switch(i) {
				case 0 : 
				case 1 : 
				case 2 : 
						 {
							test->tracking = new InterfaceTracking_Levelset3D(res, res, res, h);							

							InterfaceIntegrator_Levelset* itgLv = new InterfaceIntegrator_Levelset();

							if(i == 0)		itgLv->set_Advection(new CFDAdvection_SemiLagrangian<CFDTracking_RK2>());	
							else if(i == 1)	itgLv->set_Advection(new CFDAdvection_MacCormack(need_limiting, false));
							else if(i == 2)	itgLv->set_Advection(new CFDAdvection_BFECC(need_limiting));

							test->integrator = itgLv;

							break;
						 }

				case 3 : {					 
							int max_lv = 0; grid_index_type tmp = res;
							while(true) { if(tmp == 1) break; tmp /= 2; max_lv++; }

							test->tracking = new InterfaceTracking_SLC3D(0, res*h, max_lv);							

							test->integrator = new InterfaceIntegrator_SLC();
							break;
						 }

				case 4 : {
							InterfaceTracking_SRL3D* srl3D = new InterfaceTracking_SRL3D();
							srl3D->set(res, res, res, h, narrowband_srl, false);							
							test->tracking = srl3D;

							InterfaceIntegrator_SRL* itgSRL = new InterfaceIntegrator_SRL();
							itgSRL->set_Advection(new CFDAdvection_SemiLagrangian_SRL<CFDTracking_RK2>(need_limiting));

							test->integrator = itgSRL;
							break;
						 }
			}
			test->integrator->set_TimeStepRestriction(new CFDTimeStepRestriction(1.0e+6));

			switch(j) {
				case 0 : 
					//test->tracking->setPhi(ZalesakDistanceField3D_Enright<real>(L));
					test->tracking->setPhi(ZalesakDistanceField3D<real>(L));
					break;

				case 1 : 
					GridDistanceField3D<real> Bunny; Bunny.load("../Data/bunny_phi"); 
					Bunny.set_scale(L*0.45);
					Bunny.set_rigid_transformation(arithmetic::translate(Vector3<real>(0.5*L,0.75*L,0.5*L)));
					test->tracking->setPhi(Bunny);
					break;
			}

			test->vel = new ZalesakVelocityField3D<real>(endTime, L);
			
			while(test->frame < test->endFrame) {
				if(test->frame % 30 == 0 || (test->frame % (n_step / 4)) == 0) {
					if(i == 4)	test->tracking->createDrawableInterface(Grid3b());
					else		test->tracking->createDrawableInterface();

					char str[128];
					sprintf_s(str, 128, "%04d", test->frame);
					test->tracking->exportOBJ(output_folder + "OBJs/" + distancefield + "_" + scheme + "_" + str + ".obj");
				}				

				test->updateOneStep();				
			}

			delete test;
		}		
	}

	//test->saveImages = true;

	//test->tracking->createDrawableInterface(Grid3b());
	//scene->push_SceneObject(test);
}

void SRL_Compare_100x100_ZalesakTest(Scene* scene) {

#ifdef USE_OPENMP	
	int maxthreads = omp_get_max_threads();
	int numThreads = maxthreads;

	omp_set_num_threads(numThreads);
	log_info << endl << "OpenMP " << numThreads << " Threads are initiated" << endl << endl;
#endif

	unsigned int res = 64;
	unsigned int n_step = 628;

	const unsigned int narrowband_srl = 4;
	const unsigned int n_rotation = 1;

	T endTime = 1.0, dt; 
	T L = 1.0;

	InterfaceTrackingTest* test = new InterfaceTrackingTest();
	test->image_filename = "SRL_Zalesak100";

	log_info.set_filename("InterfaceTracking _ Zalesak _ SRL64");	test->test_name = "ZalesakSRL64";	

	dt = real(endTime) / real(n_step);
	unsigned int endframe = n_step * n_rotation + 1;
	real h = L / real(res);

	log_info << "Number of Steps _ " << n_step << endl;
	log_info << "Time Step       _ " << dt << endl;

	////////////////////////////////////////////////////////////////////////
	// Geometry & Velocity Field
	ZalesakDistanceField3D_Enright<real> zalesak(L);
	test->vel = new ZalesakVelocityField3D<real>(endTime, L);

	//////////////////////////////////////////////////////////////////////////
	// Frame
	test->FluidSimulator::set(dt, 1.0e+6);
	test->FluidSimulator::setEndFrame(endframe);	

	test->saveImages = test->saveBinaries = false;

	InterfaceTracking_SRL3D* srl3D = new InterfaceTracking_SRL3D();
	srl3D->set(res, res, res, h, narrowband_srl, false);
	srl3D->setPhi(zalesak);
	test->tracking = srl3D;

	InterfaceIntegrator_SRL* itgSRL = new InterfaceIntegrator_SRL();
	itgSRL->set_Advection(new CFDAdvection_SemiLagrangian_SRL<CFDTracking_RK2>(false));
	itgSRL->set_TimeStepRestriction(new CFDTimeStepRestriction(1.0e+6));

	test->integrator = itgSRL;	

	if(test->saveImages == false) {
		while(test->frame < test->endFrame) {

			if(test->frame % 30 == 0 || (test->frame % (n_step / 4)) == 0) {
				test->tracking->createDrawableInterface(Grid3b());

				char str[128];
				sprintf_s(str, 128, "%04d", test->frame);
				test->tracking->exportOBJ(output_folder + "OBJs/" + "SRL_Zalesak64" + "_" + str + ".obj");
			}			

			test->updateOneStep();

		}
	}

	scene->push_SceneObject(test);
}
