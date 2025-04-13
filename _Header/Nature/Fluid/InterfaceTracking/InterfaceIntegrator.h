#pragma once

#include "../FluidScene/FluidScene_Objects.h"
#include "../CFDsolver/CFDTimeStepRestriction.h"
#include "InterfaceTracking.h"

namespace fluidspace {

class NATUREDLL InterfaceIntegrator {
public : 
	InterfaceIntegrator() { numUpdates = 0; timeStepRestriction = NULL; internalBdry = NULL; }
	virtual ~InterfaceIntegrator();

	virtual void post_process() const { log_info << "InterfaceIntegrator _ Have no Post Process" << std::endl; }

	virtual void initialize(InterfaceTracking&) const = 0;

	virtual void update_InternalBoundary(const FluidScene_Objects*, const InterfaceTracking&, const CFDsolver&) 
	{ throwError("Not Implemented _ InterfaceIntegrator::update_InternalBoundary(const FluidScene_Objects*, const CFDsolver&)"); }

	virtual void update(InterfaceTracking&, CFDsolver&, real dt) = 0;
	virtual void update(InterfaceTracking&, VelocityField<real>&, real dt) = 0;

	virtual void importBinary(std::ifstream& ifs) {}
	virtual void exportBinary(std::ofstream& ofs) const {}

	void set_TimeStepRestriction(CFDTimeStepRestriction* func);

public : 
	unsigned int				numUpdates;

	CFDTimeStepRestriction*		timeStepRestriction;

	InternalBoundary*			internalBdry;
};
	
}		// fluidspace