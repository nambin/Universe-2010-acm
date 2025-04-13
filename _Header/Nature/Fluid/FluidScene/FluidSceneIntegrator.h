#ifndef _FLUIDSPACE_FLUIDSCENE_INTEGRATOR_H_
#define _FLUIDSPACE_FLUIDSCENE_INTEGRATOR_H_

#include "../../Nature.h"
#include "../Environment/Environments.h"
#include "../../../Basis/TBBProgramming/TBBProgramming.h"
#include "FluidScene_Objects.h"

namespace fluidspace {

class NATUREDLL FluidSceneIntegrator {
public : 
	FluidSceneIntegrator() {}
	FluidSceneIntegrator(real _dt0) { dt0 = _dt0; }

	//////////////////////////////////////////////////////////////////////////
	// Update
	virtual void update(real dt, FluidScene_Objects*, FluidSimulator*) const = 0;	

	virtual void post_process() const { log_info << "Fluid Scene Integrator _ Have no Post Process" << std::endl; }

public : 
	real dt0;
};

class NATUREDLL FluidSceneIntegrator_FluidOnly : public FluidSceneIntegrator {
public : 
	FluidSceneIntegrator_FluidOnly() {}
	FluidSceneIntegrator_FluidOnly(real _dt0) : FluidSceneIntegrator(_dt0) {}

	//////////////////////////////////////////////////////////////////////////
	// Update
	void update(real dt, FluidScene_Objects*, FluidSimulator*) const;	
};

class NATUREDLL FluidSceneIntegrator_StaticRigidObject : public FluidSceneIntegrator_FluidOnly {
public : 
	FluidSceneIntegrator_StaticRigidObject() : first_update(true) {}
	FluidSceneIntegrator_StaticRigidObject(real _dt0) : FluidSceneIntegrator_FluidOnly(_dt0), first_update(true) {}

	//////////////////////////////////////////////////////////////////////////
	// Update
	void update(real dt, FluidScene_Objects*, FluidSimulator*) const;	

	mutable bool first_update;
};

};			// fluidspace

#endif