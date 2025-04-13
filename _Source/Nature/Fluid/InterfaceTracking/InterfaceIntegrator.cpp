#include <stdafx.h>
#include <Nature/Fluid/InterfaceTracking/InterfaceIntegrator.h>

using namespace fluidspace;

InterfaceIntegrator::~InterfaceIntegrator() {
	if(timeStepRestriction) delete timeStepRestriction;
	 if(internalBdry) delete internalBdry;
}

void InterfaceIntegrator::set_TimeStepRestriction(CFDTimeStepRestriction* func) {
	if(timeStepRestriction) delete timeStepRestriction;
	timeStepRestriction = func;

	timeStepRestriction->print_log(log_info);
}
