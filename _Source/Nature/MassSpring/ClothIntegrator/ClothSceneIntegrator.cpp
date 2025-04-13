#include <stdafx.h>
#include <Nature/MassSpring/ClothScene/ClothSceneIntegrator.h>

using namespace springspace;

ClothSceneIntegrator::ClothSceneIntegrator() : dynamics(NULL), postDynamics(NULL), clothObjectCollision(NULL) { set(1.0f/30.0f); }
ClothSceneIntegrator::ClothSceneIntegrator(real _dt0) : dynamics(NULL), postDynamics(NULL), clothObjectCollision(NULL) { set(_dt0); }

void ClothSceneIntegrator::set(real _dt0) {
	dt0 = _dt0;
}

void ClothSceneIntegrator::init_ProvotLimiting(real _min, real _max, size_t _maxIter, bool stretch, bool shear, bool bending) {
	if(postDynamics)
		delete postDynamics;

	ProvotLimiting *_dynamics = new ProvotLimiting(stretch, shear, bending);
	_dynamics->setParameters_StrainLimiting(_min, _max, _maxIter);
	postDynamics = _dynamics;
}

void ClothSceneIntegrator::init_BridsonLimiting(real _min, real _max, size_t maxIter, real srlimit, size_t maxIter2, bool stretch, bool shear, bool bending) {
	if(postDynamics)
		delete postDynamics;

	BridsonLimiting *_dynamics = new BridsonLimiting(stretch, shear, bending);
	_dynamics->setParameters_StrainLimiting(_min, _max, maxIter);
	_dynamics->setParameters_StrainRateLimiting(srlimit, maxIter2);
	postDynamics = _dynamics;
}

void ClothSceneIntegrator::init_Bridson03_ObjectCollision(real relax) {
	if(clothObjectCollision)
		delete clothObjectCollision;

	Bridson03_Collision* collision = new Bridson03_Collision(relax);	
	clothObjectCollision = collision;
}