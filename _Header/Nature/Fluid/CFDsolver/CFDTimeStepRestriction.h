#ifndef _FLUIDSPACE_CFD_TIMESTEP_RESTRICTION_H_
#define _FLUIDSPACE_CFD_TIMESTEP_RESTRICTION_H_

#include "../../Nature.h"
#include "../Environment/Environments.h"
#include "../CFDsolver/CFDsolver2D.h"
#include "../CFDsolver/CFDsolver3D.h"

namespace fluidspace {

class NATUREDLL CFDTimeStepRestriction {
public : 
	CFDTimeStepRestriction(real _CFL) : maxCFL(_CFL) {}

	virtual void post_process() const { log_info << "CFDTimeStepRestriction _ Have no Post Process" << std::endl; }

	virtual const real		getAllowedTimeStep(const CFDsolver&, real h) const;
	virtual const real		getAllowedTimeStep(const VelocityField<real>&, real h) const;

	virtual const size_t	getNumberOfSubStep(const CFDsolver&, real dt, real h) const;
	virtual const size_t	getNumberOfSubStep(const VelocityField<real>&, real dt, real h) const;

	virtual void print_log(Log& log_) const { log_ << XmlElementName().c_str() << " _ CFL : " << maxCFL << endl; }
	static const std::string XmlElementName() { return std::string("CFDTimeStepRestriction"); }

	real maxCFL;
};

class NATUREDLL CFDTimeStepRestriction_SurfaceTension : public CFDTimeStepRestriction {
public : 
	CFDTimeStepRestriction_SurfaceTension(real _CFL, real _g, real _sigma, real _rho_neg = 75, real _rho_pos = 1) 
	: CFDTimeStepRestriction(_CFL), g(_g), sigma(_sigma), rho_neg(_rho_neg), rho_pos(_rho_pos)
	{ //if(maxCFL > 1) throwError("Error in CFDTimeStepRestriction_SurfaceTension()"); 
	}

	const real		getAllowedTimeStep(const CFDsolver&, real h) const;
	const real		getAllowedTimeStep(const VelocityField<real>&, real h) const { throwError("Not Implemented"); return 0; }

	const size_t	getNumberOfSubStep(const CFDsolver&, real dt, real h) const;
	const size_t	getNumberOfSubStep(const VelocityField<real>&, real dt, real h) const { throwError("Not Implemented"); return 0; }

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << " _ CFL : " << maxCFL << endl; }
	static const std::string XmlElementName() { return std::string("CFDTimeStepRestriction_SurfaceTension"); }
	
	real g, sigma, rho_neg, rho_pos;
};

};		// fluidspace

#endif