#ifndef _MASS_SPRING_SPACE_INTEGRATOR_H_
#define _MASS_SPRING_SPACE_INTEGRATOR_H_

#include "../../Nature.h"
#include "../Environment/Environments.h"
#include "../Core/QuadDynamicMesh.h"
#include "../Core/TriDynamicMesh.h"

namespace springspace {

class NATUREDLL MassSpringIntegrator {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	MassSpringIntegrator();
	MassSpringIntegrator(real _cfl);
	virtual ~MassSpringIntegrator() {}

	void initSimulation();

	//////////////////////////////////////////////////////////////////////////
	// Setter
	void set(real _cfl);

	//////////////////////////////////////////////////////////////////////////
	// Update	
	virtual void update(real dt, QuadDynamicMesh<real>&) = 0;
	virtual void update(real dt, TriDynamicMesh<real>&) = 0;
	virtual const real getCFL(real dt, const QuadDynamicMesh<real>&) const = 0;
	virtual const real getCFL(real dt, const TriDynamicMesh<real>&) const = 0;

public : 	
	real maxCFL;
	real g;
};

};		// namespace springspace

#endif
