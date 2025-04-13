#include <stdafx.h>
#include <Nature/Fluid/CFDsolver/CFDIntegrator.h>

using namespace std;
using namespace geometryspace;
using namespace fluidspace;

CFDIntegrator::CFDIntegrator(const bool _use_variable_dt) : use_variable_dt(_use_variable_dt) {
	addforce = NULL;
	advect = NULL;
	project = NULL;
	velBdry = NULL;
	timeStepRestriction = NULL;

	objects = NULL;

	numTimeSteps = 0;
}

CFDIntegrator::~CFDIntegrator() {
	if(addforce) delete addforce;
	if(advect) delete advect;
	if(project) delete project;
	if(velBdry) delete velBdry;
	if(timeStepRestriction) delete timeStepRestriction;
}

void CFDIntegrator::setComputingOption(const ParallelComputable& rhs) {
//	if(addforce) addforce->setComputingOption(rhs);
//	if(advect) advect->setComputingOption(rhs);
	if(project) project->setComputingOption(rhs);
//	if(velBdry) velBdry->setComputingOption(rhs);
//	if(timeStepRestriction) timeStepRestriction->setComputingOption(rhs);
}

void CFDIntegrator::set_CFDAddForce(CFDAddForce* func) {
	if(addforce) delete addforce;
	addforce = func;
//	func->setComputingOption(*this);

	addforce->print_log(log_info);
}

void CFDIntegrator::set_CFDAdvection(CFDAdvection* func) {
	if(advect) delete advect;
	advect = func;
//	func->setComputingOption(*this);

	advect->print_log(log_info);
}

void CFDIntegrator::set_CFDProjection(CFDProjection* func) {
	if(project) delete project;
	project = func;
	func->setComputingOption(*this);

	project->print_log(log_info);
}

void CFDIntegrator::set_CFDVelocityBoundaryCondition(GridBdry* u, GridBdry* v, GridBdry* w) {
	if(velBdry) delete velBdry;
	velBdry = new CFDVelocityBoundaryCondition();
	velBdry->uBdry = u; 
	velBdry->vBdry = v; 
	velBdry->wBdry = w; 

//	velBdry->setComputingOption(*this);

	velBdry->print_log(log_info);
}

void CFDIntegrator::set_CFDTimeStepRestriction(CFDTimeStepRestriction* func) {
	if(timeStepRestriction) delete timeStepRestriction;
	timeStepRestriction = func;
//	func->setComputingOption(*this);

	timeStepRestriction->print_log(log_info);
}

void CFDIntegrator::importBinary(std::ifstream& ifs) {
	if(project) project->importBinary(ifs);
}

void CFDIntegrator::exportBinary(std::ofstream& ofs) const {
	if(project) project->exportBinary(ofs);
}

void CFDIntegrator::update_substepping(CFDsolver& cfd, real dt) {
	Timer timer;

	const size_t nSubStep = timeStepRestriction->getNumberOfSubStep(cfd, dt, cfd.get_h());
	log_info << "CFD Solver SubStep : " << nSubStep << endl;

	timer.logTime("Get Number of Substep", log_info);	timer.reset();

	for(size_t i=0;i<nSubStep;++i)		
		update_one_time_step(cfd, dt / real(nSubStep));

	numTimeSteps += nSubStep;

	timer.logTime("Time for One Frame Integration for CFD Solver _ SubStepping", log_info);	
}

void CFDIntegrator::update_variable_timestep(CFDsolver& cfd, real dt) {
	Timer timer_all, timer;
	
	real dtSum = real(0.0);
	bool isEnd = false;

	while(!isEnd) {
		real dt_needed = dt - dtSum;
		real dt_allowed = timeStepRestriction->getAllowedTimeStep(cfd, cfd.get_h());

		if(dt_allowed < dt_needed)	dt_needed = dt_allowed;
		else						isEnd = true;	

		update_one_time_step(cfd, dt_needed);

		numTimeSteps++;

		dtSum += dt_needed;

		log_info << "Update Variable Time Step _ dt _ (dt_sum / dt_all) : " << dt_needed << " _ (" << dtSum << "/" << dt << ")" << endl;
		timer.printTime("Update Variable Time Step"); timer.reset();
	}

	timer_all.logTime("Time for One Frame Integration for CFD Solver _ Variable Time Step", log_info);	
}

void CFDIntegrator::update(CFDsolver& cfd, real dt) { 
	if(use_variable_dt)		update_variable_timestep(cfd, dt);
	else					update_substepping(cfd, dt); 
}

void CFDIntegrator::post_process() const {
	if(addforce) addforce->post_process();
	if(advect) advect->post_process();;
	if(project) project->post_process();;
	if(velBdry) velBdry->post_process();;
	if(timeStepRestriction) timeStepRestriction->post_process();	

	log_info << "CFDIntegrator _ Number of Updates for TimeStep _ " << numTimeSteps << std::endl;
}