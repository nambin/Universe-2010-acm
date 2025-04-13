#ifndef _MASS_SPRING_SPACE_IMPLICIT_EULER_INTEGRATOR_H_
#define _MASS_SPRING_SPACE_IMPLICIT_EULER_INTEGRATOR_H_

#include "../../../Basis/Timer/Timer.h"
#include "../MassSpringIntegrator/MassSpringIntegrator.h"
#include "../Core/MassSpringSystemMatrix.h"
#include "../Core/SpringLibrary.h"
#include "../Core/Constraints.h"

namespace springspace {

class NATUREDLL ImplicitEuler_Integrator : public MassSpringIntegrator {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	ImplicitEuler_Integrator();
	ImplicitEuler_Integrator(real _cfl);
	~ImplicitEuler_Integrator();

	void init();

	const real getCFL(real dt, const QuadDynamicMesh<real>&) const;
	const real getCFL(real dt, const TriDynamicMesh<real>&) const;

	//////////////////////////////////////////////////////////////////////////
	// Update
	void update(real dt, QuadDynamicMesh<real>&);	
	void update(real dt, TriDynamicMesh<real>&);	

	void getForce(real dt, QuadDynamicMesh<real>&);
	void getForce(real dt, TriDynamicMesh<real>&);

	//////////////////////////////////////////////////////////////////////////
	// Implicit Integration
	void allocateSystem(const QuadDynamicMesh<real>&);	
	void makeSystem(real dt, const QuadDynamicMesh<real>&);
	void makeSystem_Matrix(real dt, const QuadDynamicMesh<real>&);
	void makeSystem_RHS(real dt, const QuadDynamicMesh<real>&);

	void allocateSystem(const TriDynamicMesh<real>&);	
	void makeSystem(real dt, const TriDynamicMesh<real>&);
	void makeSystem_Matrix(real dt, const TriDynamicMesh<real>&);
	void makeSystem_RHS(real dt, const TriDynamicMesh<real>&);

	void solveConstrainedPCG(const Constraints<real>&);

public : 
	MassSpringSystemMatrix<real>		mat;
	Vec3Array<real>						b, dv, f;

	Array<std::pair<size_t, size_t>>	stretch_matIdx;
	Array<std::pair<size_t, size_t>>	shear_matIdx;
	Array<std::pair<size_t, size_t>>	bending_matIdx;
};

};		// namespace springspace

#endif
