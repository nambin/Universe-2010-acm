#include <stdafx.h>
#include <Nature/MassSpring/MassSpringIntegrator/ImplicitEuler_Integrator.h>

using namespace std;
using namespace springspace;

ImplicitEuler_Integrator::ImplicitEuler_Integrator() { init(); }
ImplicitEuler_Integrator::ImplicitEuler_Integrator(real _cfl) 
: MassSpringIntegrator(_cfl) { init(); }

ImplicitEuler_Integrator::~ImplicitEuler_Integrator() {}

void ImplicitEuler_Integrator::init() {
	MassSpringIntegrator::initSimulation();	
	mat.clear();
	b.clear();
	dv.clear();
	f.clear();
	stretch_matIdx.clear();
	shear_matIdx.clear();
	bending_matIdx.clear();
}

const real ImplicitEuler_Integrator::getCFL(real dt, const QuadDynamicMesh<real>&) const { return maxCFL; }
const real ImplicitEuler_Integrator::getCFL(real dt, const TriDynamicMesh<real>&) const  { return maxCFL; }

void ImplicitEuler_Integrator::update(real dt, QuadDynamicMesh<real>& mesh) {
	swap(mesh.pos, mesh.pos0);
	swap(mesh.vel, mesh.vel0);

	getForce(dt, mesh);
	makeSystem(dt, mesh);

	solveConstrainedPCG(mesh.constraint);
	
	for(size_t i=0;i<mesh.vel.size();++i) mesh.vel[i] = mesh.vel0[i] + dv[i];
	for(size_t i=0;i<mesh.pos.size();++i) mesh.pos[i] = mesh.pos0[i] + mesh.vel[i] * dt;	
}

void ImplicitEuler_Integrator::update(real dt, TriDynamicMesh<real>& mesh) {
	swap(mesh.pos, mesh.pos0);
	swap(mesh.vel, mesh.vel0);

	getForce(dt, mesh);
	makeSystem(dt, mesh);

	solveConstrainedPCG(mesh.constraint);

	for(size_t i=0;i<mesh.vel.size();++i) mesh.vel[i] = mesh.vel0[i] + dv[i];
	for(size_t i=0;i<mesh.pos.size();++i) mesh.pos[i] = mesh.pos0[i] + mesh.vel[i] * dt;	
}
