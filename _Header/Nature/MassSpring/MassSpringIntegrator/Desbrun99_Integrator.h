#ifndef _MASS_SPRING_SPACE_DESBRUN99_INTEGRATOR_H_
#define _MASS_SPRING_SPACE_DESBRUN99_INTEGRATOR_H_

#include "../MassSpringIntegrator/MassSpringIntegrator.h"
#include "../Core/DynamicMesh.h"
#include "../../../Principle/Arithmetic/Matrix/MatrixN.h"

namespace springspace {

class NATUREDLL Desbrun99_Integrator : public MassSpringIntegrator {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Desbrun99_Integrator();
	Desbrun99_Integrator(real _cfl);
	~Desbrun99_Integrator();

	void init();

	const real getCFL(real dt, const QuadDynamicMesh<real>&) const;
	const real getCFL(real dt, const TriDynamicMesh<real>&) const;

	//////////////////////////////////////////////////////////////////////////
	// Pre-process
	void makeFilter(real dt, const QuadDynamicMesh<real>&);
	void makeFilter(arithmetic::MatrixN<real>&, const Array<Spring<real>>&, real k) const;
	bool isFilterValid() const;

	//////////////////////////////////////////////////////////////////////////
	// Update (QuadDynamicMesh)
	void update(real dt, QuadDynamicMesh<real>&);	
	void update(real dt, TriDynamicMesh<real>&);	

	void addSpringForce(Vec3Array<real>& force, real dt, const QuadDynamicMesh<real>&) const;
	void addSpringForce(Vec3Array<real>& force, real dt, const Array<Spring<real>>&, real k, const Vec3Array<real>&, const Vec3Array<real>&) const;
	void filterForce(Vec3Array<real>& force) const;

	const Vector3<real> getAdditionalGlobalTorque(const Vec3Array<real>& force, const Vec3Array<real> &pos) const;
	void updateDynamicStates(QuadDynamicMesh<real>&, real dt, const Vec3Array<real>& force) const;
	void postCorrectionAngularMomentum(QuadDynamicMesh<real>&, real dt, const Vector3<real>& torque, const Vector3<real>& com) const;

	void postStep_IK(QuadDynamicMesh<real>&, real dt) const;

	void updateFinalVelocity(QuadDynamicMesh<real>&, real dt) const;

public : 	
	arithmetic::MatrixN<real> filter;
};

};		// namespace springspace

#endif
