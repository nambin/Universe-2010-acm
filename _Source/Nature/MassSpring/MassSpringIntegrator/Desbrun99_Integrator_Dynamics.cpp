#include <stdafx.h>
#include <vector>
#include <Nature/MassSpring/Core/Spring.h>
#include <Nature/MassSpring/Core/SpringLibrary.h>
#include <Nature/MassSpring/MassSpringIntegrator/Desbrun99_Integrator.h>
#include <Nature/MassSpring/PostProcessingDynamics/PostProcessingDynamics.h>

using namespace std;
using namespace springspace;

void Desbrun99_Integrator::addSpringForce(Vec3Array<real>& force, real dt, const QuadDynamicMesh<real>& mesh) const {
	if(mesh.stretchSpring.size() != 0)	addSpringForce(force, dt, mesh.stretchSpring, mesh.stretchFunc->k, mesh.pos0, mesh.vel0);
	if(mesh.shearSpring.size() != 0)	addSpringForce(force, dt, mesh.shearSpring, mesh.shearFunc->k, mesh.pos0, mesh.vel0);
}

void Desbrun99_Integrator::addSpringForce(Vec3Array<real>& force, real dt, const Array<Spring<real>>& springs, real k, const Vec3Array<real>& pos, const Vec3Array<real>& vel) const {
	Vector3<real> elasticforce;
	Vector3<real> dampingforce;
	real d = dt * k;

	for(size_t i=0;i<springs.size();++i) {
		const Spring<real>& s = springs[i];

		elasticforce = SpringLibrary::elasticForce(s, k, pos);
		dampingforce = SpringLibrary::dampingForce(s, d, vel);

		force[s.idx0] += elasticforce;	force[s.idx0] += dampingforce;
		force[s.idx1] -= elasticforce;	force[s.idx1] -= dampingforce;
	}
}

const Vector3<real> Desbrun99_Integrator::getAdditionalGlobalTorque(const Vec3Array<real>& force, const Vec3Array<real> &pos) const {
	Vector3<real> torque(0);
	for(size_t i=0;i<force.size();++i) 
		torque += cross(force(i), pos(i));
	return torque;
}

void Desbrun99_Integrator::updateDynamicStates(QuadDynamicMesh<real>& mesh, real dt, const Vec3Array<real>& force) const {
	Vector3<real> fext = Vector3<real>(0,g,0) * mesh.mass;
	real _mul = dt / mesh.mass;

	for(size_t i=0;i<mesh.vel.size();++i) {
		mesh.vel(i) = mesh.vel0(i) + (force(i) + fext) * _mul;
		mesh.vel(i) = mesh.constraint.getVelocity(i, mesh.vel(i));
	}

	for(size_t i=0;i<mesh.vel.size();++i) {
		mesh.pos(i) = mesh.pos0(i) + mesh.vel(i) * dt;
	}
}

void Desbrun99_Integrator::postCorrectionAngularMomentum(QuadDynamicMesh<real>& mesh, real dt, const Vector3<real>& torque, const Vector3<real>& com) const {
	real _mul = dt * dt / mesh.mass;
	for(size_t i=0;i<mesh.pos.size();++i) {
		if(!mesh.constraint.isConstrained(i)) {
			mesh.pos(i) += (cross(com - mesh.pos0[i], torque) * _mul);
		}		
	}
}

void Desbrun99_Integrator::postStep_IK(QuadDynamicMesh<real>& mesh, real dt) const {
	real ratio = 0.1f;
	size_t numIter = 200;

	ProvotLimiting provot;
	provot.setParameters_StrainLimiting(1-ratio, 1+ratio, numIter);
	provot.postProcess(dt, mesh);
}

void Desbrun99_Integrator::updateFinalVelocity(QuadDynamicMesh<real>& mesh, real dt) const {
	real dtInv = 1 / dt;
	for(size_t i=0;i<mesh.pos.size();++i) 
		mesh.vel(i) = (mesh.pos(i) - mesh.pos0(i)) * dtInv;
}