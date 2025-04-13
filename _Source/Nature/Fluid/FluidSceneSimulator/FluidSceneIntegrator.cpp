#include <stdafx.h>
#include <Nature/Fluid/FluidScene/FluidSceneIntegrator.h>

using namespace fluidspace;

void FluidSceneIntegrator_FluidOnly::update(real dt, FluidScene_Objects* objs, FluidSimulator* sim) const {
	if(sim == NULL)
		throwError("Error in FluidSceneIntegrator_FluidOnly::update()");

	if(objs != NULL)
		throwError("Error in FluidSceneIntegrator_FluidOnly::update()");
	
	real dtSum = 0;
	bool isEnd = false;

	while(!isEnd) {
		real _dt = dt - dtSum;

		if(dt0 < _dt)	_dt = dt0;
		else			isEnd = true;

		sim->update(_dt);

		dtSum += _dt;
	}	
}

void FluidSceneIntegrator_StaticRigidObject::update(real dt, FluidScene_Objects* objs, FluidSimulator* sim) const {
	if(sim == NULL)
		throwError("Error in FluidSceneIntegrator_StaticRigidObject::update()");

	if(objs == NULL) {
		FluidSceneIntegrator_FluidOnly::update(dt, NULL, sim);
		return;
	}

	if(sim->getDimension() != objs->getDimension())
		throwError("Error in FluidSceneIntegrator_StaticRigidObject::update() _ wrong dimension");

	if(first_update) {		
		sim->update_InternalBoundary(objs);
		first_update = false;
	}

	real dtSum = 0;
	bool isEnd = false;

	while(!isEnd) {
		real _dt = dt - dtSum;

		if(dt0 < _dt)	_dt = dt0;
		else			isEnd = true;

		sim->update(_dt);

		dtSum += _dt;
	}	
}
