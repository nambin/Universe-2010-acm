#pragma once

#include "../../Nature.h"

namespace fluidspace {

class FluidScene_Objects;

class NATUREDLL CFDsolver {
public : 
	CFDsolver() {}

	virtual void post_process() const { log_info << "CFDSolver _ Have no Post Process" << std::endl; }

	virtual const unsigned int getDimension() const = 0;
	virtual const real getCFL(real dt, real h) const = 0;

	virtual const real get_h() const = 0;

	virtual void update_InternalBoundary(const FluidScene_Objects&) { throwError("Not Implemented _ CFDsolver::update_InternalBoundary(const FluidScene_Objects&)"); }

	virtual const bool load(const std::string& filename) { basis::throwError("Not Implemented _ CFDsolver::load()"); return false; }
	virtual const bool save(const std::string& filename) const { basis::throwError("Not Implemented _ CFDsolver::save()"); return false; }
};

};		// fluidspace