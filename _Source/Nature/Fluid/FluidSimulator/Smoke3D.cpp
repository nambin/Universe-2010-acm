#include <stdafx.h>
#include <Nature/Fluid/FluidSimulator/Smoke3D.h>
#include <Basis/Timer/Timer.h>

using namespace fluidspace;

Smoke3D::Smoke3D() {
	log_info.set_filename("Smoke 3D");

	ioHandler = new Smoke3D_IOHandler(this);
}

Smoke3D::Smoke3D(grid_index_type nX, grid_index_type nY, grid_index_type nZ, real h, real dt, real maxCFL) {
	log_info.set_filename("Smoke 3D");

	set(nX, nY, nZ, h, dt, maxCFL);
	ioHandler = new Smoke3D_IOHandler(this);
}

Smoke3D::~Smoke3D() {
	if(integrator) 
		delete integrator;
}

void Smoke3D::setComputingOption(const ParallelComputable& rhs) {
	if(integrator)
		integrator->setComputingOption(rhs);
}

void Smoke3D::set(grid_index_type nX, grid_index_type nY, grid_index_type nZ, real h, real dt, real maxCFL) {
	FluidSimulator::set(dt, maxCFL);

	solver.set(nX, nY, nZ, h);		
	densityField.set(nX, nY, nZ, h);

	setCoefficients();

	// Setting Drawing Options
	dBoundingBox = dDensity = true;
	dVelocity = false;	
	scaleVelocity = 1.0f;	

	densityField.initialize_Shader(0,0);
}

void Smoke3D::setCoefficients(real g, real buoyancy, real rho, bool noSlip) {

	CFDIntegrator_Smoke *itg = new CFDIntegrator_Smoke();
	CFDAddForce_DensityForce* force = new CFDAddForce_DensityForce(g, buoyancy); force->setDensityTracking(&densityField);
	
	itg->set_CFDAddForce(force);
	itg->set_CFDAdvection(new CFDAdvection_SemiLagrangian<>());

	CFDProjection_SinglePhase * project = new CFDProjection_SinglePhase();
//	project->init_PoissonSolver_CG(solver); project->setPreconditioner_ICCG();
//	project->setPreconditioner_BlockRedBlackICCG(solver.nX/20,solver.nY/20,solver.nZ/20);

	project->init_PoissonSolver_GMG(solver, 10); project->set_iteration_number_GMG(3,3);
	project->set_divergence_threshold(1.0e-5);

	itg->set_CFDProjection(project);
	itg->set_CFDTimeStepRestriction(new CFDTimeStepRestriction(maxCFL));

	if(noSlip)	itg->set_CFDVelocityBoundaryCondition(new GridBoundaryCondition_U_NoSlip(), new GridBoundaryCondition_V_NoSlip(), new GridBoundaryCondition_W_NoSlip());
	else		itg->set_CFDVelocityBoundaryCondition(new GridBoundaryCondition_U_FreeSlip(), new GridBoundaryCondition_V_FreeSlip(), new GridBoundaryCondition_W_FreeSlip());

	//itg->set(&densityField, new CFDAdvection_SemiLagrangian<>());
	itg->set(&densityField, new CFDAdvection_Upwind_Conservative(1.0f));

	integrator = itg;
}

// Need to Modify
void Smoke3D::setDensity(const Grid3r &_den) {
	size_t nX = densityField.getDensity().nX, nY = densityField.getDensity().nY, nZ = densityField.getDensity().nZ;
	real dx = densityField.getDensity().dx, dy = densityField.getDensity().dy, dz = densityField.getDensity().dz, hInv = 1.0f/dx;
	real Lx = dx * nX, Ly = dy * nY, Lz = dz * nZ;

	SphereDistanceField3D<real>* sphere = new SphereDistanceField3D<real>(Vector3<real>(0.5f*Lx,0.1f*Ly,0.5f*Lz), real(0.08f*Ly));
	CubeDistanceField3D<real> *rec = new CubeDistanceField3D<real>(Vector3<real>(0.4f*Lx, 0.02f*Ly, 0.4f*Lz), Vector3<real>(0.6f*Lx, 0.04f*Ly, 0.6f*Lz));

	densityField.set(*sphere);

	((CFDIntegrator_Smoke*)integrator)->set(NULL, NULL, new DensitySource_DistanceField(sphere, 3.0));
}

void Smoke3D::update(real dt) {
	basis::Timer timer;
	log_info << endl << "frame " << frame << " at t = " << time << endl;

	integrator->update(solver, dt);

	log_info << "Density Sum : " << densityField.getTotalDensity() << endl;

	frame++;	time += dt;
	timer.logTime("Time for Simulation One Step", log_info);

	log_info << Timer::getCurrentTime().c_str() << endl;
}

void Smoke3D::update_InternalBoundary(const FluidScene_Objects* objs) {
	Timer timer;

	objects = objs;

	solver.update_InternalBoundary(*objs);
	integrator->update_InternalBoundary(objs, solver);

	timer.logTime("Smoke2D::update_InternalBoundary() Done", log_info);
}

void Smoke3D::post_process() const {
	solver.post_process();
	densityField.post_process();
	integrator->post_process();
}