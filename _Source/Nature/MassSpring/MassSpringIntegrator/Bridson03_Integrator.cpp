#include <stdafx.h>
#include <Basis/Timer/Timer.h>
#include <Nature/MassSpring/MassSpringIntegrator/Bridson03_Integrator.h>
#include <Nature/MassSpring/PostProcessingDynamics/PostProcessingDynamics.h>

using namespace std;
using namespace springspace;

Bridson03_Integrator::Bridson03_Integrator() { init(); }
Bridson03_Integrator::Bridson03_Integrator(real _cfl) 
: MassSpringIntegrator(_cfl) { init(); }

Bridson03_Integrator::~Bridson03_Integrator() {}

void Bridson03_Integrator::init() {
	MassSpringIntegrator::initSimulation();
}

const real Bridson03_Integrator::getCFL(real dt, const TriDynamicMesh<real>& mesh) const {	
	real density = mesh.getDensity();
	real _min = 1.0e+8, tmp;
	real stiff = mesh.stretchFunc->k;

	for(size_t i=0;i<mesh.stretchSpring.size();++i) {
		const Spring<real>& s = mesh.stretchSpring[i];
		tmp = length(mesh.pos[s.idx0], mesh.pos[s.idx1]) * sqrt(density/stiff);
		
		if(tmp < _min)
			_min = tmp;
	}

	return dt / _min;
}

void Bridson03_Integrator::update(real dt, TriDynamicMesh<real>& mesh) {
	swap(mesh.pos, mesh.pos0);
	swap(mesh.vel, mesh.vel0);

	velocityUpdate_Implicit(dt * 0.5f, mesh.pos0, mesh);

	// Limit Strain _ mesh.vel
	bridsonLimiting(dt, mesh);

	// Position Update
	for(size_t i=0;i<mesh.pos.size();++i) 
		mesh.pos[i] = mesh.pos0[i] + dt * mesh.vel[i];

	velocityUpdate_Explicit(dt * 0.5f, mesh.pos0, mesh);

	swap(mesh.vel, mesh.vel0);
	velocityUpdate_Implicit(dt * 0.5f, mesh.pos, mesh);

	// Limit Strain _ mesh.vel
	swap(mesh.pos, mesh.pos0);
	bridsonLimiting(dt, mesh);
	swap(mesh.pos, mesh.pos0);
}

const real Bridson03_Integrator::getCFL(real dt, const QuadDynamicMesh<real>&) const {
	basis::throwError("Not Implemented _ Bridson03_Integrator::getCFL()\n");
	return -1;
}

void Bridson03_Integrator::update(real dt, QuadDynamicMesh<real>& mesh) {
	basis::throwError("Not Implemented _ Bridson03_Integrator::update()\n");
}
