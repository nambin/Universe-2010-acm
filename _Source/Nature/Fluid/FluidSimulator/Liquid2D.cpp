#include <stdafx.h>
#include <Basis/FileIO/FileIO.h>
#include <Basis/StringLibrary/StringLibrary.h>
#include <Nature/Fluid/FluidSimulator/Liquid2D.h>
#include <Principle/Geometry/Levelset/TVD_RungeKutta.h>

using namespace fluidspace;
using namespace geometryspace;

const real narrowband_srl = 4;
const unsigned int srl_redist_freq = 5;
const bool uniform_quadrature = false;

const bool use_srl = true;
const bool use_pls = false;

Liquid2D::Liquid2D() : tracking(NULL), integrator(NULL) {
	log_info.set_filename("Liquid 2D");

	ioHandler = new Liquid2D_IOHandler(this);	
}

Liquid2D::Liquid2D(grid_index_type nX, grid_index_type nY, real h, real dt, real maxCFL, real maxCFL_Interface) : tracking(NULL), integrator(NULL) {
	log_info.set_filename("Liquid 2D");

	set(nX, nY, h, dt, maxCFL, maxCFL_Interface);
	setCoefficients();
	ioHandler = new Liquid2D_IOHandler(this);	
}

void Liquid2D::release() {	
	if(tracking)	delete tracking;	tracking = NULL;
	if(integrator)	delete integrator;	integrator = NULL;
}

Liquid2D::~Liquid2D() {
	if(tracking)	delete tracking;
	if(integrator)	delete integrator;
}

void Liquid2D::set(grid_index_type nX, grid_index_type nY, real h, real dt, real _maxCFL, real _maxCFL_Interface) {

	FluidSimulator::set(dt, _maxCFL);	
	maxCFL_Interface = _maxCFL_Interface;

	solver.set(nX, nY, h);

	if(tracking) delete tracking; tracking = NULL;

	if(use_srl == true) 
		tracking = new InterfaceTracking_SRL2D(nX, nY, h, narrowband_srl, uniform_quadrature);
	else {
		if(use_pls) tracking = new InterfaceTracking_PLS2D(nX, nY, h);	
		else		tracking = new InterfaceTracking_Levelset2D(nX, nY, h);		
	}

	// Setting Drawing Options
	dBoundingBox = dLiquid = dFilledInterface = true;
	dInternalBoundary = dGrid = dPhi = dParticle = dPressure = dVelocity = dAuxiliary = false;	
	scaleVelocity = 1.0f;	
}

void Liquid2D::setCoefficients(real _g, real _rhoPos, real _rhoNeg, bool noSlip, real narrowBand) {	
	const real sigma = 0.0728;

	assert(tracking != NULL);

	if(integrator) delete integrator;

	// CFD Solver
	CFDIntegrator_Liquid* itg = new CFDIntegrator_Liquid();
	integrator = itg;

	//itg->separation = new GridBoundaryCondition_V_CeilingSeparation(1.0);

	itg->set_CFDAddForce(new CFDAddForce_Gravity(_g));
	itg->set_CFDAdvection(new CFDAdvection_SemiLagrangian<>());
	
//	CFDProjection_Liquid* project = new CFDProjection_FreeSurface();	
//	itg->set_CFDTimeStepRestriction(new CFDTimeStepRestriction(maxCFL));	
	
	CFDProjection_Liquid* project = new CFDProjection_FreeSurface_Gibou(false, sigma);
	itg->set_CFDTimeStepRestriction(new CFDTimeStepRestriction_SurfaceTension(1.0f, _g, 0));

	project->set_InterfaceTracking(tracking);
//	project->init_PoissonSolver_CG(solver);	project->setPreconditioner_ICCG(); 
//	project->setPreconditioner_BlockRedBlackICCG(solver.nX / 20, solver.nY / 20);
	project->init_PoissonSolver_GMG(solver, 10); project->set_iteration_number_GMG(3,5);
	project->set_divergence_threshold(1.0e-3);

	itg->set_CFDProjection(project);

	if(noSlip)	itg->set_CFDVelocityBoundaryCondition(new GridBoundaryCondition_U_NoSlip(), new GridBoundaryCondition_V_NoSlip(), new GridBoundaryCondition_W_NoSlip());
	else		itg->set_CFDVelocityBoundaryCondition(new GridBoundaryCondition_U_FreeSlip(), new GridBoundaryCondition_V_FreeSlip(), new GridBoundaryCondition_W_FreeSlip());

	// Interface Integrator
	InterfaceIntegrator* itg_interface(NULL);
	if(use_srl) {
		InterfaceIntegrator_SRL* itg_srl = new InterfaceIntegrator_SRL();

		SRL_Redistancing* redist = new SRL_Redistancing_PDE_WENO(narrowband_srl, 0.5, true, 0);
		//SRL_Redistancing* redist = new SRL_Redistancing_Direct(narrowband_srl);		

		redist->set_redistance_cfl_frequency(8);
		redist->set_redistance_min_max_updates(5,50);
		//redist->set_redistance_min_max_updates(0,0);

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
//		itg_lv->set_Redistancing(new Levelset_Redistancing_PDE_ENO<real,RungeKutta_Euler>(narrowBand, 0.5f));

		itg_interface = itg_lv;
	}

	itg_interface->set_TimeStepRestriction(new CFDTimeStepRestriction(maxCFL_Interface));

	Levelset_Extrapolation<real>* extrapolation = NULL;
	extrapolation = new Levelset_Extrapolation_PDE<real>(narrowBand, 0.5f);

	itg->set(tracking, itg_interface, extrapolation);	
}

void Liquid2D::setPhi(const Grid2r &_phi) {
	const real Lx = tracking->get_Lx();
	const real Ly = tracking->get_Ly();
	const real h = tracking->get_h();

	RectangleDistanceField2D<real> damBreaking(Vector2<real>(-10.0f * Lx), Vector2<real>(0.3f*Lx, 0.7*Ly));
	tracking->setPhi(damBreaking);

//	WaterDropDistanceField2D<real> drop(Vector2<real>(0.5f*Lx, 0.75f*Ly), 0.12f*Ly, 0.1f*Ly);
//	tracking->setPhi(drop);

	//CircleDistanceField2D<real> circle(Vector2<real>(0.51f*Lx, 0.5f*Ly), 0.5f);
	//tracking->setPhi(circle);

	integrator->integrator->initialize(*tracking);
	tracking->createDrawableInterface(solver.isKinematicBoundary);	
}

void Liquid2D::setPhi(const DistanceField2D<real>& field) {
	tracking->setPhi(field);
	integrator->integrator->initialize(*tracking);	
	tracking->createDrawableInterface(solver.isKinematicBoundary); 
}

void Liquid2D::update(real _dt)	{
	Timer timer;

	log_info << endl << "Frame " << frame << " at t = " << time << endl;

	integrator->update(solver, _dt);

	frame++; time += _dt;
	timer.logTime("Time per One Step Simulation", log_info);	timer.reset();

	log_info << "Volume : " << tracking->getTotalVolume(solver.isKinematicBoundary) << endl;
	
	solver.printMaxVel(log_info);
	tracking->createDrawableInterface(solver.isKinematicBoundary);

	timer.logTime("Time for creating Drawable Interface", log_info);
	log_info << Timer::getCurrentTime().c_str() << endl;
}

void Liquid2D::loadBinary(std::size_t frame) {
	char filename[512];

	sprintf_s( filename, (output_folder + get_concat_filename("FluidBinary/Liquid2D _ Simulator ", frame)).c_str(), frame );
	std::ifstream ifs(StringLibrary::string2wstring(filename).c_str(), ios::binary);

	FluidSimulator::importBinary(ifs);
	::importBinary(ifs, maxCFL_Interface);

	integrator->importBinary(ifs);

	sprintf_s( filename, (output_folder + get_concat_filename("FluidBinary/Liquid2D _ " + solver.name_class(), frame)).c_str(), frame );
	solver.load(filename);

	sprintf_s( filename, (output_folder + get_concat_filename("FluidBinary/Liquid2D _ " + tracking->name_class(), frame)).c_str(), frame );
	tracking->load(filename);

	tracking->createDrawableInterface(solver.isKinematicBoundary);
}

void Liquid2D::saveBinary(std::size_t frame) const {
	char filename[512];

	sprintf_s( filename, (output_folder + get_concat_filename("FluidBinary/Liquid2D _ Simulator ", frame)).c_str(), frame );
	std::ofstream ofs(StringLibrary::string2wstring(filename).c_str(), ios::binary);

	FluidSimulator::exportBinary(ofs);
	::exportBinary(ofs, maxCFL_Interface);

	integrator->exportBinary(ofs);
	
	sprintf_s( filename, (output_folder + get_concat_filename("FluidBinary/Liquid2D _ " + solver.name_class(), frame)).c_str(), frame );
	solver.save(filename);

	sprintf_s( filename, (output_folder + get_concat_filename("FluidBinary/Liquid2D _ " + tracking->name_class(), frame)).c_str(), frame );
	tracking->save(filename);
}

void Liquid2D::saveRenderObject(std::size_t frame) const {
	char filename[512];

	sprintf_s( filename, (output_folder + get_concat_filename("FluidBinary/Liquid2D _ " + tracking->name_class(), frame)).c_str(), frame );
	tracking->save_compressed(filename);
}

void Liquid2D::update_InternalBoundary(const FluidScene_Objects* objs) {
	Timer timer;

	objects = objs;
	solver.update_InternalBoundary(*objs);
	integrator->update_InternalBoundary(objs, solver);

	timer.logTime("Liquid2D::update_InternalBoundary() Done", log_info);
}

void Liquid2D::post_process() const {
	solver.post_process();
	tracking->post_process();
	integrator->post_process();
}