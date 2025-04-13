#ifndef _NATURE_MASS_SPRING_CLOTH_INTEGRATOR_H_
#define _NATURE_MASS_SPRING_CLOTH_INTEGRATOR_H_

#include "../../Nature.h"
#include "../Environment/Environments.h"
#include "../MassSpringIntegrator/Desbrun99_Integrator.h"
#include "../MassSpringIntegrator/ImplicitEuler_Integrator.h"
#include "../MassSpringIntegrator/Bridson03_Integrator.h"
#include "../PostProcessingDynamics/PostProcessingDynamics.h"
#include "../Collision/ClothObjectCollision.h"

namespace springspace {

class NATUREDLL ClothSceneIntegrator {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	ClothSceneIntegrator();
	ClothSceneIntegrator(real dt0);
	virtual ~ClothSceneIntegrator() {}

	//////////////////////////////////////////////////////////////////////////
	// Setter
	void set(real dt0);

	template<typename MESH>
	void init_ImplicitEuler_Integrator(const MESH&, real cfl);

	template<typename MESH>
	void init_Bridson03_Integrator(const MESH&, real cfl);

	void init_ProvotLimiting(real min, real max, size_t iter, bool stretch = true, bool shear = true, bool bending = false);
	void init_BridsonLimiting(real min, real max, size_t iter, real srlimit, size_t iter2, bool stretch = true, bool shear = true, bool bending = false);

	void init_Bridson03_ObjectCollision(real relax);

	//////////////////////////////////////////////////////////////////////////
	// Update
	virtual void update(real dt, ClothScene_Objects*, DynamicMesh*) = 0;	

public : 
	real	dt0;

	MassSpringIntegrator*		dynamics;
	PostProcessingDynamics*		postDynamics;
	ClothObjectCollision*		clothObjectCollision;
};

template<typename MESH>
void ClothSceneIntegrator::init_ImplicitEuler_Integrator(const MESH& mesh, real maxCFL) {
	if(mesh.size() == 0)
		basis::throwError("Should initialize DynamicMesh first\n");

	if(dynamics)
		delete dynamics;
	
	ImplicitEuler_Integrator* _dynamics = new ImplicitEuler_Integrator(maxCFL);
	_dynamics->allocateSystem(mesh);
	dynamics = _dynamics;
}

template<typename MESH>
void ClothSceneIntegrator::init_Bridson03_Integrator(const MESH& mesh, real maxCFL) {
	if(mesh.size() == 0)
		basis::throwError("Should initialize DynamicMesh first\n");

	if(dynamics)
		delete dynamics;

	Bridson03_Integrator* _dynamics = new Bridson03_Integrator(maxCFL);
	_dynamics->allocateSystem(mesh);
	dynamics = _dynamics;
}

class NATUREDLL ClothSceneIntegrator_Basic : public ClothSceneIntegrator {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	ClothSceneIntegrator_Basic(real dt);

	//////////////////////////////////////////////////////////////////////////
	// Update
	virtual void update(real dt, ClothScene_Objects*, DynamicMesh*);
	
protected : 
	void update(real dt, ClothScene_Objects*, TriDynamicMesh<real>*);
	void update(real dt, ClothScene_Objects*, QuadDynamicMesh<real>*);
};

};		// springspace

#endif