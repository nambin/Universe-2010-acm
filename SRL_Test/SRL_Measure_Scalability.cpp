#include <stdafx.h>
#include "SRL_Experiment.h"

using namespace std;
using namespace geometryspace;
using namespace rigidspace;
using namespace fluidspace;

typedef fluidspace::real T;

void SRL_Measure_Scalability(Scene* scene) {

	unsigned int res, endframe;

	std::cout << "Input the Domain Resolution (100~200) : ";	std::cin >> res;
	std::cout << "Input the Frame Number (5~20) : ";			std::cin >> endframe;

	int maxthreads = 1;

#ifdef USE_OPENMP	
	maxthreads = omp_get_max_threads();
#endif	

	const unsigned int narrowband_srl = 5;
	const unsigned int n_rotation = 1;

	unsigned int n_step = 628;
	T endTime = 100, dt = real(endTime) / real(n_step);
	T L = 100;
	T h = L / real(res);	

	InterfaceTrackingTest* test = new InterfaceTrackingTest();
	test->saveImages = test->saveBinaries = false;

	for(int i=1;i<=maxthreads;++i) {
	//for(int i=maxthreads;i>=1;--i) {

		if(i != 1 && i != 4 && i != 8) 			
		//if(i != 1 && i != 2 && i != 4 && i != 7 && i != 8) 		
			continue;

		char str[256];
		sprintf_s(str, 256, "%d", i);
		log_info.set_filename(string("MeasureScalability_") + str);
		test->test_name = string("Parallel") + str;

		test->release();

#ifdef USE_OPENMP	
		omp_set_num_threads(i);
#endif

		log_info << endl << "OpenMP " << i << " Threads are initiated" << endl << endl;
		
		////////////////////////////////////////////////////////////////////////
		// Geometry & Velocity Field
		GridDistanceField3D<real> Bunny; Bunny.load("../Data/bunny_phi"); 
		Bunny.set_scale(L*0.7);
		Bunny.set_rigid_transformation(arithmetic::translate(Vector3<real>(0.5*L,0.5*L,0.5*L)));

		test->vel = new ZalesakVelocityField3D<real>(endTime, L);

		//////////////////////////////////////////////////////////////////////////
		// Frame
		test->FluidSimulator::set(dt, 1.0e+6);
		test->FluidSimulator::setEndFrame(endframe);	

		InterfaceTracking_SRL3D* srl3D = new InterfaceTracking_SRL3D();
		srl3D->set(res, res, res, h, narrowband_srl, false);
		srl3D->setPhi(Bunny);
		test->tracking = srl3D;

		InterfaceIntegrator_SRL* itgSRL = new InterfaceIntegrator_SRL();
		itgSRL->set_Advection(new CFDAdvection_SemiLagrangian_SRL<CFDTracking_RK4>(true));
		itgSRL->set_Redistancing(new SRL_Redistancing_PDE_WENO(narrowband_srl, 0.5, true, 0));
		itgSRL->set_TimeStepRestriction(new CFDTimeStepRestriction(1.0e+6));

		test->integrator = itgSRL;

		test->tracking->createDrawableInterface(Grid3b());

		while(test->frame < test->endFrame) 
			test->updateOneStep();

		test->post_process();
	}

	scene->push_SceneObject(test);
}