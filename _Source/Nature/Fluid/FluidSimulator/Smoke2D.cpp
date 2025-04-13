#include <stdafx.h>
#include <Basis/Timer/Timer.h>
#include <Nature/Fluid/FluidScene/FluidScene_Objects.h>
#include <Nature/Fluid/FluidSimulator/Smoke2D.h>
#include <Nature/Fluid/CFDsolver/CFDIntegrator.h>

using namespace fluidspace;

Smoke2D::Smoke2D() {
	log_info.set_filename("Smoke 2D");
	
	ioHandler = new Smoke2D_IOHandler(this);
}

Smoke2D::Smoke2D(grid_index_type nX, grid_index_type nY, real h, real dt, real maxCFL) {
	log_info.set_filename("Smoke 2D");

	set(nX, nY, h, dt, maxCFL);
	ioHandler = new Smoke2D_IOHandler(this);
}

Smoke2D::~Smoke2D() {
	if(integrator) 
		delete integrator;
}

void Smoke2D::setComputingOption(const ParallelComputable& rhs) {
	if(integrator)
		integrator->setComputingOption(rhs);
}

void Smoke2D::set(grid_index_type nX, grid_index_type nY, real h, real dt, real maxCFL)	{
	FluidSimulator::set(dt, maxCFL);

	solver.set(nX, nY, h);	
	densityField.set(nX, nY, h);	

	setCoefficients();

	// Setting Drawing Options
	dBoundingBox = dDensity = true;
	dKinematicBoundary = dPressure = dVelocity = false;	
	scaleVelocity = 1.0f;
}

void Smoke2D::setCoefficients(real g, real buoyancy, real rho, bool noSlip) {
	
	CFDIntegrator_Smoke *itg = new CFDIntegrator_Smoke();
	CFDAddForce_DensityForce* force = new CFDAddForce_DensityForce(g, buoyancy); force->setDensityTracking(&densityField);

	itg->set_CFDAddForce(force);
	itg->set_CFDAdvection(new CFDAdvection_SemiLagrangian<>());

	CFDProjection_SinglePhase * project = new CFDProjection_SinglePhase();
//	project->init_PoissonSolver_CG(solver);	//project->setPreconditioner_ICCG();
//	project->setPreconditioner_BlockRedBlackICCG(solver.nX / 20, solver.nY / 20);	

	project->init_PoissonSolver_GMG(solver, 10); project->set_iteration_number_GMG(2,2);
	project->set_divergence_threshold(1.0e-3);

	itg->set_CFDProjection(project);
	itg->set_CFDTimeStepRestriction(new CFDTimeStepRestriction(maxCFL));

	if(noSlip)	itg->set_CFDVelocityBoundaryCondition(new GridBoundaryCondition_U_NoSlip(), new GridBoundaryCondition_V_NoSlip(), new GridBoundaryCondition_W_NoSlip());
	else		itg->set_CFDVelocityBoundaryCondition(new GridBoundaryCondition_U_FreeSlip(), new GridBoundaryCondition_V_FreeSlip(), new GridBoundaryCondition_W_FreeSlip());

//	itg->set(&densityField, new CFDAdvection_SemiLagrangian<>());
	itg->set(&densityField, new CFDAdvection_Upwind_Conservative(1.0f));

	integrator = itg;
}

void Smoke2D::setDensity(const Grid2r &den) {
	const real Lx = densityField.density.dx * densityField.density.nX, Ly = densityField.density.dy * densityField.density.nY; 

	CircleDistanceField2D<real>* circle = new CircleDistanceField2D<real>(Vector2<real>(0.5f*Lx,0.04f*Ly), real(0.02f*Ly));
	RectangleDistanceField2D<real> *rec = new RectangleDistanceField2D<real>(Vector2<real>(0.4f*Lx, 0.02f*Ly), Vector2<real>(0.6f*Lx, 0.04f*Ly));

	densityField.set(*circle);

	((CFDIntegrator_Smoke*)integrator)->set(NULL, NULL, new DensitySource_DistanceField(circle, 3.0));
}

void Smoke2D::update(real dt) {
	basis::Timer timer;
	log_info << endl << "Frame " << frame << " at t = " << time << endl;

	integrator->update(solver, dt);

	log_info << "Density Sum : " << densityField.getTotalDensity() << endl;

	frame++; time += dt;
	timer.logTime("Time for Simulation One Step", log_info);

	log_info << Timer::getCurrentTime().c_str() << endl;
}

void Smoke2D::update_InternalBoundary(const FluidScene_Objects* objs) {
	Timer timer;

	objects = objs;

	solver.update_InternalBoundary(*objs);
	integrator->update_InternalBoundary(objs, solver);

	timer.logTime("Smoke2D::update_InternalBoundary() Done", log_info);
}

void Smoke2D::post_process() const {
	solver.post_process();
	densityField.post_process();
	integrator->post_process();
}