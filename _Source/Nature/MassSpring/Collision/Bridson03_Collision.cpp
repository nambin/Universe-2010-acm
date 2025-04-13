#include <stdafx.h>
#include <Nature/MassSpring/Collision/ClothObjectCollision.h>

using namespace std;
using namespace springspace;

Bridson03_Collision::Bridson03_Collision(real relax) {
	assert(relax > 0);
	relaxTime = relax;	
}

Bridson03_Collision::~Bridson03_Collision() {
}

void Bridson03_Collision::treatCollision(const ClothScene_Objects& scene, real dt, TriDynamicMesh<real>& mesh) {
	treatCollision(scene, dt, mesh.pos0, mesh.vel, mesh.constraint);
	for(size_t i=0;i<mesh.pos.size();++i) 
		mesh.pos[i] = mesh.pos0[i] + mesh.vel[i] * dt;	
}

void Bridson03_Collision::treatCollision(const ClothScene_Objects& scene, real dt, springspace::QuadDynamicMesh<real>& mesh) {
	treatCollision(scene, dt, mesh.pos0, mesh.vel, mesh.constraint);
	for(size_t i=0;i<mesh.pos.size();++i) 
		mesh.pos[i] = mesh.pos0[i] + mesh.vel[i] * dt;	
}

void Bridson03_Collision::treatCollision(const ClothScene_Objects& scene, real dt, const Vec3Array<real>& pos0, Vec3Array<real>& vel, Constraints<real>& constraints) const {

	real relaxTimeInv = 1.0f / relaxTime;

	real phi, phinew, mu;
	Vector3<real> vel_o, normal;

	real vn, vn_o, vn_new;
	Vector3<real> vel_t, vel_t_o, vel_t_rel, vel_t_rel_new;	

	ClothScene_Objects::QueryInfo query;

	real eps = 1.0e-3;
	real epsContact = eps;
	for(size_t i=0;i<pos0.size();++i) {
		scene.query(query, pos0[i]);

		phi = query.phi;
		vel_o = query.vel;
		normal = query.normal;
		mu = query.mu;

		phinew = phi + dt * dot(vel[i] - vel_o, normal);
		phinew -= eps;
		
		if(phinew > 0)
			continue;
		
		vn = dot(vel[i], normal);	vel_t = vel[i] - normal * vn;
		vn_o = dot(vel_o, normal);	vel_t_o = vel_o - normal * vn_o;

		vn_new = vn - phinew * relaxTimeInv;

		vel_t_rel = vel_t - vel_t_o;
		vel_t_rel_new = MAX(real(0), 1.0f - mu * (vn_new - vn) / mag(vel_t_rel)) * vel_t_rel;

		//vel[i] = vn_new * normal + vel_t;		
		vel[i] = vn_new * normal + vel_t_o + vel_t_rel_new;

		constraints.addConstraint(i, normal);		
		constraints.setVelocity(i, vel[i]);
	}	
}	