#pragma once

#include "../CFDsolver/CFDAdvection.h"
#include "../InterfaceTracking/InterfaceTracking_SRL2D.h"
#include "../InterfaceTracking/InterfaceTracking_SRL3D.h"
#include "../InterfaceTracking/SRL_Redistancing.h"
#include "InterfaceIntegrator.h"

namespace fluidspace {

class NATUREDLL InterfaceIntegrator_SRL : public InterfaceIntegrator {
public : 
	InterfaceIntegrator_SRL() : advect(NULL), redistancing(NULL), time_for_tracking(0) {}
	virtual ~InterfaceIntegrator_SRL();

	void post_process() const;

	void initialize(InterfaceTracking&) const {}

	void update_InternalBoundary(const FluidScene_Objects*, const InterfaceTracking&, const CFDsolver&);

	void update(InterfaceTracking&, CFDsolver&, real dt);
	void update(InterfaceTracking&, VelocityField<real>&, real dt);

	template<class SRL, class CFD> void update_srl(SRL&, CFD&, real dt);

	void set_Advection(CFDAdvection*);
	void set_Redistancing(SRL_Redistancing*);

	void importBinary(std::ifstream& ifs);
	void exportBinary(std::ofstream& ofs) const;

public : 
	CFDAdvection*			advect;
	SRL_Redistancing*		redistancing;

	real					time_for_tracking;
};

template<class SRL, class CFD> void InterfaceIntegrator_SRL::update_srl(SRL& srl, CFD& cfd, real _dt) {

	srl.srl->print_debug_info();

	const real internalbdry_extrapolate_band = timeStepRestriction->maxCFL + 2;

	basis::Timer timer_all;
	basis::Timer timer;
	
	const real C_cfl = cfd.getCFL(1.0f, srl.srl->h);
	const real cfl = _dt * C_cfl;

	log_info << "Interface CFL : " << cfl << endl;

	const size_t nSubStep = size_t(cfl/timeStepRestriction->maxCFL + 1);
	const real dt = _dt / nSubStep;
	
	log_info << "SRL SubStep : " << nSubStep << endl;
	timer.logTime("Get Number of Substep", log_info);	timer.reset();

	basis::Timer timer_tracking;

	for(size_t i=0;i<nSubStep;i++) {
		advect->advect(srl, cfd, dt);
		timer.logTime("Advect SRL", log_info);	timer.reset();

		if(internalBdry) {
			internalBdry->extrapolate_phi(srl, internalbdry_extrapolate_band);	
			timer.logTime("Extrapolate SRL", log_info); timer.reset();
		}	

		if(redistancing && redistancing->is_times_up(dt, C_cfl)) {
			redistancing->redistance(srl);	
			timer.logTime("Redistacing SRL", log_info);	timer.reset();

			if(internalBdry) {
				internalBdry->extrapolate_phi(srl, internalbdry_extrapolate_band);	
				timer.logTime("Extrapolate SRL", log_info); timer.reset();
			}	
		}			
	}	

//	srl.srl0->phi = srl.srl->phi;		// synchronize coarse grid phi
										// for safety

	srl.srl->print_debug_info();
	timer_all.logTime("Total Time for SRL Integration", log_info);

	time_for_tracking += timer_tracking.getTime();
}

};			// fluidspace