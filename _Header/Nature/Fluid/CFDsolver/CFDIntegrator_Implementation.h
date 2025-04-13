#pragma once

namespace fluidspace {

template<class CFD> 
void CFDIntegrator_Smoke::update_one_time_step_template(CFD& cfd, real dt) {
	Timer timer;

	if(source)	{ 
		source->add(field, dt); 
		timer.logTime("Add Density", log_info);	timer.reset();
	}

	if(addforce) {
		addforce->addForce(cfd, dt);
		velBdry->enforce(cfd);
		timer.logTime("Add Force", log_info);	timer.reset();
	}

	if(project) {
		project->project(cfd, dt);
		velBdry->enforce(cfd);
		timer.logTime("Project", log_info);	timer.reset();
	}	

	if(field->getDimension() == 2)  { advect_density->advect((DensityTracking2D&)(*field), (CFDsolver2D&)cfd, dt); }
	else							{ advect_density->advect((DensityTracking3D&)(*field), (CFDsolver3D&)cfd, dt); }

	timer.logTime("Advect Density", log_info);	timer.reset();

	advect->advect(cfd, dt, *velBdry);
	timer.logTime("Advect Velocity", log_info);	timer.reset();
}

template<class CFD> 
void CFDIntegrator_Liquid::update_one_time_step_template(CFD& cfd, real dt) {
	Timer timer;
	
	if(advect) {
		advect->advect(cfd, dt, *velBdry);
		timer.logTime("Advect Velocity", log_info);	timer.reset();
	}

	if(addforce) {
		addforce->addForce(cfd, dt);	
		velBdry->enforce(cfd);
		timer.logTime("Add Force", log_info);	timer.reset();
	}		

	cfd.setVelocity_InternalBoundary(objects);
	velBdry->enforce(cfd);
	timer.logTime("Velocity Boundary Condition _ Internal & External", log_info);	timer.reset();	

	if(source)	{ 
		source->add(tracking, cfd, dt); 		
		timer.logTime("Add Liquid", log_info);	timer.reset();
	}

	if(project) {
		project->project(cfd, dt);
		velBdry->enforce(cfd);
		timer.logTime("Project", log_info);	timer.reset();			
	}		

	extrapolateVelocity(cfd);
	timer.logTime("Velocity Extrapolation", log_info);	timer.reset();

	cfd.setVelocity_InternalBoundary(objects);
	velBdry->enforce(cfd);
	timer.logTime("Velocity Boundary Condition _ Internal & External", log_info);	timer.reset();	

	if(separation)
		separation->enforce(cfd.v);

	integrator->update(*tracking, cfd, dt);
	timer.logTime("Interface Update", log_info);	timer.reset();	

	velBdry->enforce(cfd);
}

};		// fluidspace