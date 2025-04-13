#pragma once

namespace fluidspace {

template<class LV, class CFD> void InterfaceIntegrator_Levelset::update_lv(LV& levelset, CFD& cfd, real _dt) {

	const real band = timeStepRestriction->maxCFL + 2;

	basis::Timer timer;

	const size_t nSubStep = timeStepRestriction->getNumberOfSubStep(cfd, _dt, levelset.h);
	const real dt = _dt / nSubStep;
		
	log_info << "Levelset SubStep : " << nSubStep << endl;
	timer.logTime("Get Number of Substep", log_info);	timer.reset();

	for(size_t i=0;i<nSubStep;i++) {
		advect->advect(levelset, cfd, dt);				
		timer.logTime("Advect Levelset", log_info);	timer.reset();

		if(internalBdry) {
			internalBdry->extrapolate_phi(levelset,band);	
			timer.logTime("Extrapolate Levelset", log_info); timer.reset();
		}		
	}	

	if(redistancing) {
		redistancing->redistance(levelset);				
		timer.logTime("Redistacing Levelset", log_info); timer.reset();
		
		if(internalBdry) {
			internalBdry->extrapolate_phi(levelset,band);	
			timer.logTime("Extrapolate Levelset", log_info); timer.reset();
		}		
	}	
}

template<class PLS, class CFD> void InterfaceIntegrator_PLS::update_pls(PLS& pls, CFD& cfd, real _dt) {
	
	basis::Timer timer, timer_s;
	
	const size_t nSubStep = timeStepRestriction->getNumberOfSubStep(cfd, _dt, pls.h);	
	const real dt = _dt / nSubStep;
	
	log_info << "ParticleLevelset SubStep : " << nSubStep << endl;
	timer.logTime("Get Number of Substep", log_info);	timer.reset();

	basis::Timer timer_all;

	double adv_lv = 0, adv_particle = 0, error_crt = 0, redist = 0;
	for(size_t i=0;i<nSubStep;i++) {		

		timer_s.reset();
		
		advect->advect(pls, cfd, dt);
		adv_lv += timer_s.getTime(); timer_s.reset();

		particleAdvect->advect(pls, cfd, dt);
		adv_particle += timer_s.getTime(); timer_s.reset();

		errCorrector->errorCorrect(pls, true);
		error_crt += timer_s.getTime(); timer_s.reset();

		particleMnp->removeDeepParticles(pls);	
	}

	timer.logTime("Advect PLS", log_info);	timer.reset();

	if(redistancing) {
		timer_s.reset();

		redistancing->redistance(pls);
		redist += timer_s.getTime(); timer_s.reset();

		errCorrector->errorCorrect(pls, true);
		error_crt += timer_s.getTime(); timer_s.reset();
	}	

	timer.logTime("Redistacing PLS", log_info);	timer.reset();

	particleMnp->removeDeepParticles(pls);	

	particleMnp->adjustRadius(pls);

	timer_all.logTime("Total Time for PLS", log_info);

	if(numUpdates % particleMnp->reseed_cycle == 0)
		particleMnp->reseed(pls);

	timer.logTime("Post-Processing PLS", log_info);	timer.reset();

	log_info << "Time for Advecting Levelset phi : " << adv_lv << endl;
	log_info << "Time for Advecting Particle     : " << adv_particle << endl;
	log_info << "Time for Correcting Error       : " << error_crt << endl;
	log_info << "Time for Redistancing phi       : " << redist << endl;

	pls.printParticleInformation(log_info);
}

};			// fluidspace