#include <stdafx.h>
#include <Basis/Random/RandomGen.h>
#include "SRL_Experiment.h"

using namespace std;
using namespace geometryspace;
using namespace rigidspace;
using namespace fluidspace;

typedef fluidspace::real T;

class LiquidSource_Falling_Armadillo : public LiquidSource {
public : 
	LiquidSource_Falling_Armadillo() {
		armadillo.load("../Data/armadillo_phi"); 
	}		

	void set_vertical_falling(real domain_size_xz, real domain_size_y, real ts, real te, real scale_min, real scale_max, real v, int num_per_sec) {
		const int num = (te-ts) * num_per_sec;

		t0_position.resize(num);
		t0_scale.resize(num);
		t0_time_point.resize(num);
		vel0 = Vector3<real>(0,v,0);

		std::vector<int> range; 
		range.push_back(0);	range.push_back(3);
		range.push_back(2);	range.push_back(1);
		range.push_back(3);	range.push_back(1);
		range.push_back(0);	range.push_back(2);

		for(size_t i=8;i<num;++i)
			range.push_back(range[i%8]);

		for(unsigned int i=0;i<num;++i) {
			const real t_interval = (te - ts) / real(num);

			if(i==0) t0_time_point[i] = ts + 0.001;
			else	 t0_time_point[i] = (random.getFloat() * t_interval) + (i * t_interval) + ts;

			t0_scale[i] = (random.getFloat() * (scale_max - scale_min) + scale_min) * domain_size_xz;

			const real domain_4 = domain_size_xz * 0.25;

			t0_position[i][0] = domain_4 + (range[i]%2) * domain_4 * 2.0;
			t0_position[i][2] = domain_4 + (range[i]/2) * domain_4 * 2.0;
			t0_position[i][1] = domain_size_y - t0_scale[i] * 0.6;
		}
	}

	void set_across_falling(real domain_size_xz, real domain_size_y, real ts, real te, real scale_min, real scale_max, real vx, real vy) {

		const unsigned int num = 4;
		const unsigned int num_collide = 3;

		t1_position.resize(num+num_collide*2);
		t1_scale.resize(num+num_collide*2);
		t1_time_point.resize(num+num_collide*2);	
		t1_rotation.resize(num+num_collide*2);	

		//vel1 = Vector3<real>(vx,vy,0);
		const real ts_2_interval = (te-ts) * real(num*0.5) / real(num*0.5 + num_collide);

		const real domain_4 = domain_size_xz * 0.25;
		const real domain_3 = domain_size_xz * 0.333;

		std::vector<int> range; 
		range.push_back(2);	range.push_back(1);		
		range.push_back(3);	range.push_back(1);		
		range.push_back(2);	range.push_back(3);		

		real t_interval = (ts_2_interval) / real(num);
		for(unsigned int i=0;i<num;++i) {		
			t1_time_point[i] = (0.5 * t_interval) + (i * t_interval) + ts;
			t1_scale[i] = (random.getFloat() * (scale_max - scale_min) + scale_min) * domain_size_xz;			

			t1_position[i][0] = domain_4 + (range[i]%2) * domain_4 * 2.0;
			t1_position[i][2] = t1_scale[i] * 0.5 + (range[i]%2) * domain_4 + (domain_size_xz - t1_scale[i]) * 0.5 * random.getFloat();
			t1_position[i][1] = domain_size_y - t1_scale[i] * 0.5;

			if(t1_position[i][0] > domain_size_xz * 0.5)	vel1.push_back(Vector3<real>(-vx, vy, 0));
			else											vel1.push_back(Vector3<real>(+vx, vy, 0));

			if(t1_position[i][0] > domain_size_xz * 0.5)	{ t1_rotation[i] = arithmetic::rotateZ(real(+135)) * arithmetic::rotateY(real(-90)); } //vel1[0] = - abs(vel1[0]); }
			else											{ t1_rotation[i] = arithmetic::rotateZ(real(-135)) * arithmetic::rotateY(real(+90)); } //vel1[0] = + abs(vel1[0]); }
		}

		t_interval = (te - ts - ts_2_interval) / real(num_collide);
		for(unsigned int i=0;i<num_collide;++i) {

			if(i == 0) {
				t1_time_point[num+2*i] = t1_time_point[num+2*i+1] = (random.getFloat() * t_interval) + i * t_interval + ts_2_interval + ts;
				t1_time_point[num+2*i+1] += (t_interval * 0.3 * (random.getFloat() * 2.0f - 1.0f));

				t1_scale[num+2*i  ] = (random.getFloat() * (scale_max - scale_min) + scale_min) * domain_size_xz;
				t1_scale[num+2*i+1] = (random.getFloat() * (scale_max - scale_min) + scale_min) * domain_size_xz;

				t1_position[num+2*i  ][0] = domain_4 ;
				t1_position[num+2*i+1][0] = domain_4 * 3.0 ;

				t1_position[num+2*i  ][2] = real(range[i]*2-1) * 0.5 * domain_3 ;
				t1_position[num+2*i+1][2] = real(range[i]*2-1) * 0.5 * domain_3 ;
				t1_position[num+2*i  ][1] = domain_size_y - t1_scale[num+2*i  ] * 0.6;
				t1_position[num+2*i+1][1] = domain_size_y - t1_scale[num+2*i+1] * 0.6;

				vel1.push_back(Vector3<real>(+vx*0.3, vy, 0));
				vel1.push_back(Vector3<real>(-vx*0.3, vy, 0));

				if(t1_position[num+2*i][0] > domain_size_xz * 0.5)	{ t1_rotation[num+2*i] = arithmetic::rotateZ(real(+135)) * arithmetic::rotateY(real(-90)); } //vel1[0] = - abs(vel1[0]); }
				else												{ t1_rotation[num+2*i] = arithmetic::rotateZ(real(-135)) * arithmetic::rotateY(real(+90)); } //vel1[0] = + abs(vel1[0]); }

				if(t1_position[num+2*i+1][0] > domain_size_xz * 0.5)	{ t1_rotation[num+2*i+1] = arithmetic::rotateZ(real(+135)) * arithmetic::rotateY(real(-90)); } //vel1[0] = - abs(vel1[0]); }
				else													{ t1_rotation[num+2*i+1] = arithmetic::rotateZ(real(-135)) * arithmetic::rotateY(real(+90)); } //vel1[0] = + abs(vel1[0]); }
			}
			else if(i==1) {

				t1_time_point[num+2*i] = t1_time_point[num+2*i+1] = (0.5 * t_interval) + i * t_interval + ts_2_interval + ts;
				t1_time_point[num+2*i+1] += (t_interval * 0.1 * (random.getFloat() * 2.0f - 1.0f));

				t1_scale[num+2*i  ] = (0.5 * (scale_max - scale_min) + scale_min) * domain_size_xz;
				t1_scale[num+2*i+1] = (0.5 * (scale_max - scale_min) + scale_min) * domain_size_xz;

				t1_position[num+2*i  ][0] = domain_4 ;
				t1_position[num+2*i+1][0] = domain_4 * 3.0 ;

				t1_position[num+2*i  ][2] = domain_4 * 3.0 ;
				t1_position[num+2*i+1][2] = domain_4 ;
				t1_position[num+2*i  ][1] = domain_size_y - t1_scale[num+2*i  ] * 0.6;
				t1_position[num+2*i+1][1] = domain_size_y - t1_scale[num+2*i+1] * 0.6;

				vel1.push_back(Vector3<real>(+vx*0.4, vy, -vx*0.4));
				vel1.push_back(Vector3<real>(-vx*0.4, vy, +vx*0.4));

				if(t1_position[num+2*i][0] > domain_size_xz * 0.5)	{ t1_rotation[num+2*i] = arithmetic::rotate(real(+135), unitVector(Vector3<real>(1,0,1))) * arithmetic::rotateY(real(-45)); } //vel1[0] = - abs(vel1[0]); }
				else												{ t1_rotation[num+2*i] = arithmetic::rotate(real(-135), unitVector(Vector3<real>(1,0,1))) * arithmetic::rotateY(real(+135)); } //vel1[0] = + abs(vel1[0]); }

				if(t1_position[num+2*i+1][0] > domain_size_xz * 0.5)	{ t1_rotation[num+2*i+1] = arithmetic::rotate(real(+135), unitVector(Vector3<real>(1,0,1))) * arithmetic::rotateY(real(-45)); } //vel1[0] = - abs(vel1[0]); }
				else													{ t1_rotation[num+2*i+1] = arithmetic::rotate(real(-135), unitVector(Vector3<real>(1,0,1))) * arithmetic::rotateY(real(+135)); } //vel1[0] = + abs(vel1[0]); }

			}
			else if(i==2) {

				t1_time_point[num+2*i] = t1_time_point[num+2*i+1] = (0.5 * t_interval) + i * t_interval + ts_2_interval + ts;
				t1_time_point[num+2*i+1] += (t_interval * 0.2 * (random.getFloat() * 2.0f - 1.0f));

				t1_scale[num+2*i  ] = (0.5 * (scale_max - scale_min) + scale_min) * domain_size_xz;
				t1_scale[num+2*i+1] = (0.5 * (scale_max - scale_min) + scale_min) * domain_size_xz;

				t1_position[num+2*i  ][0] = domain_4 ;
				t1_position[num+2*i+1][0] = domain_4 * 3.0 ;

				t1_position[num+2*i  ][2] = domain_4 ;
				t1_position[num+2*i+1][2] = domain_4 * 3.0 ;
				t1_position[num+2*i  ][1] = domain_size_y - t1_scale[num+2*i  ] * 0.6;
				t1_position[num+2*i+1][1] = domain_size_y - t1_scale[num+2*i+1] * 0.6;

				vel1.push_back(Vector3<real>(+vx*0.4, vy, +vx*0.4));
				vel1.push_back(Vector3<real>(-vx*0.4, vy, -vx*0.4));

				if(t1_position[num+2*i][0] > domain_size_xz * 0.5)	{ t1_rotation[num+2*i] = arithmetic::rotateY(real(180)) * arithmetic::rotate(real(+135), unitVector(Vector3<real>(1,0,-1))) * arithmetic::rotateY(real(+45)); } //vel1[0] = - abs(vel1[0]); }
				else												{ t1_rotation[num+2*i] = arithmetic::rotateY(real(180)) * arithmetic::rotate(real(-135), unitVector(Vector3<real>(1,0,-1))) * arithmetic::rotateY(real(-135)); } //vel1[0] = + abs(vel1[0]); }

				if(t1_position[num+2*i+1][0] > domain_size_xz * 0.5) { t1_rotation[num+2*i+1] = arithmetic::rotateY(real(180)) * arithmetic::rotate(real(+135), unitVector(Vector3<real>(1,0,-1))) * arithmetic::rotateY(real(+45)); } //vel1[0] = - abs(vel1[0]); }
				else												 { t1_rotation[num+2*i+1] = arithmetic::rotateY(real(180)) * arithmetic::rotate(real(-135), unitVector(Vector3<real>(1,0,-1))) * arithmetic::rotateY(real(-135)); } //vel1[0] = + abs(vel1[0]); }

			}
		}		
	}

	void add(InterfaceTracking* tracking, CFDsolver& cfd, real dt) const {

		for(size_t i=0;i<t0_time_point.size();i++) {
			if(time <= t0_time_point[i] && t0_time_point[i] < time + dt) {								
				armadillo.set_rigid_transformation(arithmetic::translate(t0_position[i]) * arithmetic::rotateZ(real(180.0)));
				armadillo.set_scale(t0_scale[i]);
				tracking->addSource(armadillo);

				CFDsolver3D& cfd3D = (CFDsolver3D&)cfd;
				setVel(cfd3D, vel0, armadillo);		
			}
		}

		Transform<real> rotation;
		for(size_t i=0;i<t1_time_point.size();++i) {
			if(time <= t1_time_point[i] && t1_time_point[i] < time + dt) {		

				CFDsolver3D& cfd3D = (CFDsolver3D&)cfd;

				//if(t1_position[i][0] > cfd3D.Lx * 0.5)	{ rotation = arithmetic::rotateZ(real(+135)) * arithmetic::rotateY(real(-90)); } //vel1[0] = - abs(vel1[0]); }
				//else									{ rotation = arithmetic::rotateZ(real(-135)) * arithmetic::rotateY(real(+90)); } //vel1[0] = + abs(vel1[0]); }

				armadillo.set_rigid_transformation(arithmetic::translate(t1_position[i]) * t1_rotation[i]);
				armadillo.set_scale(t1_scale[i]);
				tracking->addSource(armadillo);
				
				setVel(cfd3D, vel1[i], armadillo);		
			}
		}
	
		time += dt;
	}	
	
protected : 
	basis::RandomGen<real>				random;

	std::vector<Vector3<real>>			t0_position;
	std::vector<real>					t0_scale;
	std::vector<real>					t0_time_point;
	Vector3<real>						vel0;

	std::vector<Vector3<real>>			t1_position;
	std::vector<real>					t1_scale;
	std::vector<real>					t1_time_point;
	std::vector<Transform<T>>			t1_rotation;
	std::vector<Vector3<real>>			vel1;
	
	//mutable Vector3<real>				vel1;

	mutable GridDistanceField3D<real>	armadillo;
};

void SRL_FluidSimulation_3D_Demo_02(Scene* scene) {

	const real g = - 9.81, sigma = 0.000;

	const int endFrame = 270;
	const real dt = 1.0 / 30.0;	
	const real cfl = 3.5, cfl_interface = 1.5;

	const real narrowband_srl = 4;
	const unsigned int redist_freq = 15;
	const unsigned int redist_min = 5, redist_max = 50;

	const real ratio = 0.5;
	const int nX = 48 * ratio, nY = 64 * ratio, nZ = 48 * ratio;
	const real h = 2.5 / double(nX);	

	const real Lx = nX*h, Ly = nY*h, Lz = nZ*h;

	FluidScene *fluidscene = new FluidScene(dt);
	fluidscene->initFluidSceneIntegrator(new FluidSceneIntegrator_StaticRigidObject(dt));	
	fluidscene->setEndFrame(endFrame);
	fluidscene->enableSaveBinary(500);
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

	LiquidSource_Falling_Armadillo* source = new LiquidSource_Falling_Armadillo();
//	source->set_vertical_falling(Lx, Ly, 0.1, 3, 0.5, 0.5, -4.0);
//	source->set_across_falling(Lx, Ly, 3.1, 6, 0.5, 0.5, 3.0, -3.0);

	source->set_vertical_falling(Lx, Ly, 0.01, 5.5, 0.5, 0.5, -4.0, 5);
	source->set_across_falling(Lx, Ly, 5.6, 7.5, 0.5, 0.5, 3.0, -3.0);
	
	itg->set(liquid3D->tracking, itg_srl, extrapolation, source);	
	liquid3D->integrator = itg;

	////////////////////////////////////////////////////////////////////////	
	// Set Phi
	JointDistanceField3D<real> jointfield;

	WaterDropDistanceField3D<real> water(Vector3<real>(10000.0), 0.0, 0.025*Ly);
	jointfield.add_distancefields(&water);

	liquid3D->setPhi(jointfield);

	fluidscene->initFluidSimulator(liquid3D);	

	////////////////////////////////////////////////////////////////////////	
	// Set Object
	//FluidScene_Objects3D* object3D = new FluidScene_Objects3D();
	//object3D->push_DynamicObject(StaticRigidBody3D<T>(Sphere<T>(Vector3<T>(nX*h*0.5f, nY*h*0.3f, nZ*h*0.5f), nY*h*0.15f)));
	//fluidscene->initFluidScene_Objects(object3D);
}