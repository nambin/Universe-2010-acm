#ifndef _FLUIDSPACE_LEVELSET_INTEGRATOR_H_
#define _FLUIDSPACE_LEVELSET_INTEGRATOR_H_

#include "../../../Principle/Geometry/Levelset/Levelset_Redistancing.h"
#include "../../../Principle/Geometry/Levelset/Levelset_Extrapolation.h"
#include "../CFDsolver/CFDAdvection.h"
#include "../InterfaceTracking/InterfaceTracking_Levelset2D.h"
#include "../InterfaceTracking/InterfaceTracking_Levelset3D.h"
#include "../InterfaceTracking/InterfaceTracking_PLS2D.h"
#include "../InterfaceTracking/ParticleLevelset_Mnp.h"
#include "../InterfaceTracking/ParticleLevelset_ParticleAdvect.h"
#include "../InterfaceTracking/ParticleLevelset_ErrorCorrection.h"
#include "../FluidLibrary/InternalBoundary.h"
#include "InterfaceIntegrator.h"

namespace fluidspace {

using namespace geometryspace;

class NATUREDLL InterfaceIntegrator_Levelset : public InterfaceIntegrator {
public : 
	InterfaceIntegrator_Levelset();
	virtual ~InterfaceIntegrator_Levelset();

	void initialize(InterfaceTracking&) const {}

	void set_Advection(CFDAdvection*);
	void set_Redistancing(Levelset_Redistancing<real>*);

	void update_InternalBoundary(const FluidScene_Objects*, const InterfaceTracking&, const CFDsolver&);

	void update(InterfaceTracking&, CFDsolver&, real dt);
	void update(InterfaceTracking&, VelocityField<real>&, real dt);

	template<class LV, class CFD> void update_lv(LV&, CFD&, real dt);

	virtual void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; }
	
	static const std::string XmlElementName() { return std::string("InterfaceIntegrator_Levelset"); }

public : 
	CFDAdvection*					advect;
	Levelset_Redistancing<real>*	redistancing;	
};

class NATUREDLL InterfaceIntegrator_PLS : public InterfaceIntegrator_Levelset {
public : 
	InterfaceIntegrator_PLS();
	~InterfaceIntegrator_PLS();

	void initialize(InterfaceTracking&) const;
	
	void update(InterfaceTracking&, CFDsolver&, real dt);
	void update(InterfaceTracking&, VelocityField<real>&, real dt);

	template<class PLS, class CFD> void update_pls(PLS&, CFD&, real dt);

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; }

	static const std::string XmlElementName() { return std::string("InterfaceIntegrator_PLS"); }

public : 
	PLS_ParticleMnp*		particleMnp;
	PLS_PartlcleAdvector*	particleAdvect;
	PLS_ErrorCorrection*	errCorrector;
};

};		// fluidspace

#include "InterfaceIntegrator_Levelset_Integration.h"

#endif