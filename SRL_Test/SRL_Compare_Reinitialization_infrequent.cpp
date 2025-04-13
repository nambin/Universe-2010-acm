#include <stdafx.h>
#include "SRL_Experiment.h"

using namespace std;
using namespace geometryspace;
using namespace rigidspace;
using namespace fluidspace;

typedef fluidspace::real T;

void SRL_Compare_Reinitialization_RotatingBunny_Infrequent(Scene* scene) {

#ifdef USE_OPENMP	
	int maxthreads = omp_get_max_threads();
	int numThreads = maxthreads;

	omp_set_num_threads(numThreads);
	log_info << endl << "OpenMP " << numThreads << " Threads are initiated" << endl << endl;
#endif

	const bool clipping = true;
	const real rotate_angle = 0;

	const unsigned int redist_freq = 15;
	const unsigned int redist_min = 5, redist_max = 50;	
	
	const unsigned int narrowband_srl = 5;
	const unsigned int n_rotation = 20;
	const unsigned int res = 128;
	const unsigned int endTime = 200; 

	const grid_index_type nX(res), nY(res), nZ(res);
	const real L = 100;
	const real h = L / real(res);

//	VelocityField3D<T>* vel_field = new RotationVelocityField3D<real>(endTime, L, Vector3<T>(0.5*L, 0.5*L, 0.5*L), arithmetic::rotateY(rotate_angle).transformVector(Vector3<T>(0,0,1)));
	VelocityField3D<T>* vel_field = new ZalesakVelocityField3D<real>(endTime, L);
	
	CFDTimeStepRestriction time_restrict(1.0f);
	size_t n_step = time_restrict.getNumberOfSubStep(*vel_field, endTime, h) - 1;

	const T dt = T(endTime) / (T(n_step));

	log_info << "Number of Steps for one rotation : " << n_step << endl << endl;

	const unsigned int endframe = int(endTime / dt + 1) * n_rotation + n_rotation + 1;	

	GridDistanceField3D<real> bunny; bunny.load("../Data/bunny_phi"); 
	Transform<T> transform;
	//Transform<T> transform =  arithmetic::translate(+ Vector3<real>(0.5*L,0,0.5*L)) 
	//						* arithmetic::rotateY(rotate_angle)
	//						* arithmetic::translate(- Vector3<real>(0.5*L,0,0.5*L)) ;

	bunny.set_scale(L*0.45);	
	bunny.set_rigid_transformation(transform * arithmetic::translate(Vector3<real>(0.5*L,0.75*L,0.5*L)));
	
	std::string scheme;
	InterfaceTrackingTest test;

	int i = -1;
	cout << "0 : WENO3rd" << endl;
	cout << "1 : WENO" << endl;
	cout << "2 : Direct" << endl;
	cout << "3 : None" << endl;
	cout << endl;

	while( i <= -1 || 4 <= i ) {
		cout << "Select : ";
		cin >> i;
	}
	
	//for(int i=0;i<1;++i) 
	{
		test.release();

		if(i==0)	scheme = "WENO3rd";		
		if(i==1)	scheme = "WENO";
		if(i==2)	scheme = "Direct";
		if(i==3)	scheme = "None";

		log_info.set_filename("InterfaceTracking_Rotating Bunny_" + scheme);	test.test_name = "RotatingBunny" + scheme;

		////////////////////////////////////////////////////////////////////////
		// Geometry & Velocity Field			
		test.vel = vel_field->clone();

		//////////////////////////////////////////////////////////////////////////
		// Frame
		test.FluidSimulator::set(dt, 1.0e+6);
		test.FluidSimulator::setEndFrame(endframe);	
		test.saveBinaries = true;
		test.saveImages = false;		

		InterfaceTracking_SRL3D* srl3D = new InterfaceTracking_SRL3D();
		srl3D->set(nX, nY, nZ, h, narrowband_srl, false);
		srl3D->setPhi(bunny);
		test.tracking = srl3D;

		for(int k=0;k<5;++k)
			cout << srl3D->srl->r_sub[k] << endl;		

		InterfaceIntegrator_SRL* itgSRL = new InterfaceIntegrator_SRL();
		itgSRL->set_Advection(new CFDAdvection_SemiLagrangian_SRL<CFDTracking_RK4>(clipping));

		SRL_Redistancing* redist = NULL;

		if( i == 0 )		redist = new SRL_Redistancing_PDE_WENO3rd(narrowband_srl, 0.4, true, 0);		
		else if( i == 1 )	redist = new SRL_Redistancing_PDE_WENO(narrowband_srl, 0.4, true, 0);	
		else if( i == 2 )	redist = new SRL_Redistancing_Direct(narrowband_srl);
		else if( i == 3 )	redist = NULL;

		if(redist != NULL) {
			redist->set_redistance_cfl_frequency(redist_freq);	
			redist->set_redistance_min_max_updates(redist_min, redist_max);	
		}		

		if(redist == NULL) {
			redist = new SRL_Redistancing_PDE_WENO(narrowband_srl, 0.4, true, 0);	
			redist->set_redistance_cfl_frequency(1.0e+6);	
			redist->set_redistance_min_max_updates(5, unsigned int(n_step));	
		}

		itgSRL->set_Redistancing(redist);		

		itgSRL->set_TimeStepRestriction(new CFDTimeStepRestriction(1.0e+6));
		test.integrator = itgSRL;

		while(test.frame < test.endFrame) {			
			if((test.frame % n_step == 0) || (test.frame % (n_step / 3)) == 0) {
				test.tracking->createDrawableInterface(Grid3b());

				char str[128];
				sprintf_s(str, 128, "%04d", test.frame);
				test.tracking->exportOBJ(output_folder + "OBJs/" + "Bunny" + "_" + scheme + "_" + str + ".obj", transform.getInverse());
			}		

			//if(test.frame % n_step == 0) {
			//	InterfaceTracking_SRL3D* srl = dynamic_cast<InterfaceTracking_SRL3D*>(test.tracking);
			//	SRL_Redistancing_PDE_WENO(narrowband_srl, 0.4, true, 0).redistance(*srl);
			//}

			test.updateOneStep();
		}		
	}	
}
