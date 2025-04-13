#include <stdafx.h>
#include <Nature/MassSpring/MassSpringIntegrator/ImplicitEuler_Integrator.h>

using namespace std;
using namespace springspace;

void ImplicitEuler_Integrator::getForce(real dt, QuadDynamicMesh<real>& mesh) {
	Vector3<real> fext(0,g,0);
	for(size_t i=0;i<f.size();++i)
		f[i] = fext * mesh.mass;

	if(mesh.stretchSpring.size() != 0)	SpringLibrary::addSpringForce(f, mesh.stretchSpring, mesh.stretchFunc, mesh.pos0, mesh.vel0);
	if(mesh.shearSpring.size() != 0)	SpringLibrary::addSpringForce(f, mesh.shearSpring, mesh.shearFunc, mesh.pos0, mesh.vel0);
	if(mesh.bendingSpring.size() != 0)	SpringLibrary::addSpringForce(f, mesh.bendingSpring, mesh.bendingFunc, mesh.pos0, mesh.vel0);
}

void ImplicitEuler_Integrator::getForce(real dt, TriDynamicMesh<real>& mesh) {
	Vector3<real> fext(0,g,0);
	for(size_t i=0;i<f.size();++i)
		f[i] = fext * mesh.mass[i];

	if(mesh.stretchSpring.size() != 0)	SpringLibrary::addSpringForce(f, mesh.stretchSpring, mesh.stretchFunc, mesh.pos0, mesh.vel0);
}