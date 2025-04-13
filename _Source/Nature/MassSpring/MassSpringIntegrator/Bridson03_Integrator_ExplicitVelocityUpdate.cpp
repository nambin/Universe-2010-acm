#include <stdafx.h>
#include <Nature/MassSpring/MassSpringIntegrator/Bridson03_Integrator.h>

using namespace std;
using namespace springspace;

void Bridson03_Integrator::velocityUpdate_Explicit(real dt, const Vec3Array<real> &x, TriDynamicMesh<real>& mesh) {
	Vector3<real> fext(0,g,0);
	for(size_t i=0;i<f.size();++i)
		f[i] = fext;

	if(mesh.stretchSpring.size() != 0)	SpringLibrary::addSpringForce(f, mesh.stretchSpring, mesh.stretchFunc, x, mesh.vel0);

	mesh.constraint.filterVelocities(f);
	for(size_t i=0;i<mesh.vel.size();++i) 
		mesh.vel[i] = mesh.vel0[i] + dt * f[i];	
}