#ifndef _MASS_SPRING_SPACE_BRIDSON03_INTEGRATOR_H_
#define _MASS_SPRING_SPACE_BRIDSON03_INTEGRATOR_H_

#include "../MassSpringIntegrator/MassSpringIntegrator.h"
#include "../PostProcessingDynamics/PostProcessingDynamics.h"

namespace springspace {

class NATUREDLL Bridson03_Integrator : public MassSpringIntegrator {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Bridson03_Integrator();
	Bridson03_Integrator(real _cfl);
	~Bridson03_Integrator();

	void init();

	const real getCFL(real dt, const QuadDynamicMesh<real>&) const;
	const real getCFL(real dt, const TriDynamicMesh<real>&) const;

	//////////////////////////////////////////////////////////////////////////
	// Update
	void update(real dt, QuadDynamicMesh<real>&);	
	void update(real dt, TriDynamicMesh<real>&);	

	void velocityUpdate_Implicit(real dt, const Vec3Array<real> &x, TriDynamicMesh<real>&);	
	void velocityUpdate_Explicit(real dt, const Vec3Array<real> &x, TriDynamicMesh<real>&);
	void bridsonLimiting(real dt, TriDynamicMesh<real>&);

	void velocityUpdate_Implicit_(real dt, const Vec3Array<real> &x, TriDynamicMesh<real>&);

	//////////////////////////////////////////////////////////////////////////
	// System
	void allocateSystem(const TriDynamicMesh<real>& mesh);

	void makeSystem_RHS(real dt, const TriDynamicMesh<real>&);

	bool pcgSolve(Vec3Array<real>& ret, const Vec3Array<real>& b, const TriDynamicMesh<real>&, real dt, real &resi, size_t &iter);
	void makeDiagonalPreconditioner(const TriDynamicMesh<real>&, real dt);
	void applyDiagonalPreconditioner(Vec3Array<real>& ret, const Vec3Array<real>& b) const;	
	void mul(Vec3Array<real> &ret, const TriDynamicMesh<real>&, real dt, const Vec3Array<real> &v);
	
	void makeSystem(real dt, const TriDynamicMesh<real>&);
	void makeSystem_Matrix(real dt, const TriDynamicMesh<real>&);

	void updateSystem_Matrix(MassSpringSystemMatrix<real>& mat, const Vec3Array<real> &pos0, const Vec3Array<real> &vel0, real dt, 
							const Array<Spring<real>>& spring, const SpringFunctor<real>* func, const Array<std::pair<size_t, size_t>>& matIdx);

public : 
	MassSpringSystemMatrix<real>		mat;

	Vec3Array<real>						b, dv, f;
	Vec3Array<real>						z, r, s;
	Array<real>							invDiag;
	real								tolerance;
	size_t								maxIter;

	Array<std::pair<size_t, size_t>>	stretch_matIdx;
};

};		// springspace

#endif