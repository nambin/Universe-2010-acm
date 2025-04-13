#include <stdafx.h>
#include <Nature/Fluid/FluidSimulator/Liquid3D.h>

using namespace fluidspace;
using namespace geometryspace;

const real narrowband_srl = 4;
const unsigned int srl_redist_freq = 8;
const bool uniform_quadrature = false;

const bool use_srl = true;
const bool use_pls = false;

Liquid3D::Liquid3D() {
	log_info.set_filename("Liquid 3D");

	tracking = NULL;
	integrator = NULL; 

	ioHandler = new Liquid3D_IOHandler(this);	
}

Liquid3D::Liquid3D(grid_index_type nX, grid_index_type nY, grid_index_type nZ, real h, real dt, real maxCFL, real maxCFL_Interface) : tracking(NULL), integrator(NULL) {
	log_info.set_filename("Liquid 3D");

	set(nX, nY, nZ, h, dt, maxCFL, maxCFL_Interface);
//	setCoefficients();
	ioHandler = new Liquid3D_IOHandler(this);	
}

Liquid3D::~Liquid3D() {	
	release();
}

void Liquid3D::release() {	
	if(tracking)	delete tracking;	tracking = NULL;
	if(integrator)	delete integrator;	integrator = NULL;
}

void Liquid3D::set(grid_index_type nX, grid_index_type nY, grid_index_type nZ, real h, real dt, real _maxCFL, real _maxCFL_Interface) {

	FluidSimulator::set(dt, _maxCFL);
	maxCFL_Interface = _maxCFL_Interface;

	solver.set(nX, nY, nZ, h);

	if(tracking) delete tracking; tracking = NULL;

	if(use_srl == true) 
		tracking = new InterfaceTracking_SRL3D(nX, nY, nZ, h, narrowband_srl, uniform_quadrature);
	else {
		if(use_pls) tracking = new InterfaceTracking_PLS3D(nX, nY, nZ, h);	
		else		tracking = new InterfaceTracking_Levelset3D(nX, nY, nZ, h);		
	}
	
	// Setting Drawing Options
	dBoundingBox = dSurface = true;
	dMesh = dInnerBoundary = false;
}

void Liquid3D::setCoefficients(real _g, real _rhoPos, real _rhoNeg, bool noSlip, real narrowBand) {
//	_g = 0;
	const real sigma = 0.0005;
	const real rho_neg = 300, rho_pos = 1.226;

	assert(tracking != NULL);
	
	if(integrator) delete integrator;

	// CFD Solver
	CFDIntegrator_Liquid* itg = new CFDIntegrator_Liquid();
	integrator = itg;

	itg->set_CFDAddForce(new CFDAddForce_Gravity(_g));
//	itg->set_CFDAddForce(new CFDAddForce_SurfaceTension(_g, sigma, rho_neg, rho_pos));

	itg->set_CFDAdvection(new CFDAdvection_SemiLagrangian<>());

//	CFDProjection_Liquid* project = new CFDProjection_FreeSurface();
//	itg->set_CFDTimeStepRestriction(new CFDTimeStepRestriction(maxCFL));

	CFDProjection_Liquid* project = new CFDProjection_FreeSurface_Gibou(true, sigma);
	itg->set_CFDTimeStepRestriction(new CFDTimeStepRestriction_SurfaceTension(1.0f, _g, sigma));	

	//CFDProjection_Liquid* project = new CFDProjection_CSF(rho_neg, rho_pos);
	//itg->set_CFDTimeStepRestriction(new CFDTimeStepRestriction_SurfaceTension(10.0f, _g, sigma, rho_neg, rho_pos));

	project->set_InterfaceTracking(tracking);
//	project->init_PoissonSolver_CG(solver);	project->setPreconditioner_ICCG();
//	project->setPreconditioner_BlockRedBlackICCG(solver.nX / 20, solver.nY / 20, solver.nZ / 20);
	project->init_PoissonSolver_GMG(solver, 12); project->set_iteration_number_GMG(3,5);
	project->set_divergence_threshold(1.0e-5);

	itg->set_CFDProjection(project);	

	if(noSlip)	itg->set_CFDVelocityBoundaryCondition(new GridBoundaryCondition_U_NoSlip(), new GridBoundaryCondition_V_NoSlip(), new GridBoundaryCondition_W_NoSlip());
	else		itg->set_CFDVelocityBoundaryCondition(new GridBoundaryCondition_U_FreeSlip(), new GridBoundaryCondition_V_FreeSlip(), new GridBoundaryCondition_W_FreeSlip());

	// Interface Integrator
	InterfaceIntegrator* itg_interface(NULL);
	if(use_srl) {
		InterfaceIntegrator_SRL* itg_srl = new InterfaceIntegrator_SRL();

		SRL_Redistancing* redist = new SRL_Redistancing_PDE_WENO(narrowband_srl, 0.5, true, 0);
		//SRL_Redistancing* redist = new SRL_Redistancing_PDE_2ndUpwind(narrowband_srl, 0.5, true, 0);

		redist->set_redistance_cfl_frequency(srl_redist_freq);
		redist->set_redistance_min_max_updates(5,50);

		itg_srl->set_Advection(new CFDAdvection_SemiLagrangian_SRL<CFDTracking_RK2>(true));		
		itg_srl->set_Redistancing(redist);	

		itg_interface = itg_srl;
	}
	else {
		InterfaceIntegrator_Levelset* itg_lv(NULL);

		if(use_pls) itg_lv = new InterfaceIntegrator_PLS();
		else		itg_lv = new InterfaceIntegrator_Levelset();

		itg_lv->set_Advection(new CFDAdvection_SemiLagrangian<CFDTracking_RK2>());
		itg_lv->set_Redistancing(new Levelset_RedistancingFMM<real>(narrowBand));
//		itg_lv->set_Redistancing(new Levelset_Redistancing_PDE_WENO<real,RungeKutta_Euler>(narrowBand, 0.5f));

		itg_interface = itg_lv;
	}
	
	itg_interface->set_TimeStepRestriction(new CFDTimeStepRestriction(maxCFL_Interface));
	
	Levelset_Extrapolation<real>* extrapolation = NULL;
	if(dynamic_cast<CFDProjection_CSF*>(project) == NULL)
		extrapolation = new Levelset_Extrapolation_PDE<real>(narrowBand, 0.5f);
		 
	itg->set(tracking, itg_interface, extrapolation);	
}

void Liquid3D::setPhi(const Grid3r &_phi) {
	const real Lx = tracking->get_Lx(), Ly = tracking->get_Ly(), Lz = tracking->get_Lz(), h = tracking->get_h();

	WaterDropDistanceField3D<real> drop(Vector3<real>(0.5f*Lx, 0.83*Ly, 0.5f*Lz), 0.12f*Ly, 0.15f*Ly);
	tracking->setPhi(drop);

//	SphereDistanceField3D<real> sphere(Vector3<real>(0.5f*Lx, 0.5f*Ly, 0.5f*Lz), 1.0);
//	tracking->setPhi(sphere);
	
	integrator->integrator->initialize(*tracking);	

//	if(saveImages) {
		tracking->createDrawableInterface(solver.isKinematicBoundary); 
//	}
}

void Liquid3D::setPhi(const DistanceField3D<real>& field) {
	tracking->setPhi(field);
	integrator->integrator->initialize(*tracking);	

//	if(saveImages) {
		tracking->createDrawableInterface(solver.isKinematicBoundary); 
//	}
}

void Liquid3D::update(real _dt)	{
	Timer timer;

	log_info << endl << "Frame " << frame << " at t = " << time << endl;

	integrator->update(solver, _dt);

	frame++; time += _dt;
	timer.logTime("Time per One Step Simulation", log_info); timer.reset();

	log_info << "Volume : " << tracking->getTotalVolume(solver.isKinematicBoundary) << endl;

	timer.logTime("Time for get Volume", log_info); timer.reset();

	if(saveImages) {
		tracking->createDrawableInterface(solver.isKinematicBoundary); 
		timer.logTime("Time for creating Drawable Interface", log_info);
	}	

	log_info << Timer::getCurrentTime().c_str() << endl;
}

void Liquid3D::update_InternalBoundary(const FluidScene_Objects* objs) {
	Timer timer;

	objects = objs;
	solver.update_InternalBoundary(*objs);
	integrator->update_InternalBoundary(objs, solver);

	timer.logTime("Liquid3D::update_InternalBoundary() Done", log_info);
}

void Liquid3D::saveRenderObject(std::size_t frame) const {
	char filename[512];

	sprintf_s( filename, (output_folder + get_concat_filename("FluidBinary/Liquid3D _ " + tracking->name_class(), frame)).c_str(), frame );
	tracking->save_compressed(filename);
}

void Liquid3D::post_process() const {
	solver.post_process();
	tracking->post_process();
	integrator->post_process();
}