#include <stdafx.h>
#include <Nature/MassSpring/MassSpringIntegrator/Bridson03_Integrator.h>

using namespace std;
using namespace springspace;

void Bridson03_Integrator::bridsonLimiting(real dt, TriDynamicMesh<real>& mesh) {
	BridsonLimiting limiting;

	limiting.setParameters_StrainLimiting(0.9f, 1.1f, 15);
	limiting.setParameters_StrainRateLimiting(0.1f, 20);

	std::vector<Array<Spring<real>>*> springs;
	springs.push_back(&mesh.stretchSpring);

	real maxStr = limiting.strainLimiting(mesh.pos0, mesh.vel, dt, mesh.constraint, springs);
	real maxStrRate = limiting.strainRateLimiting(mesh.pos0, mesh.vel, dt, mesh.constraint, springs);

//	cout << "Max Strain      : " << maxStr << endl;
//	cout << "Max Strain Rate : " << maxStrRate << endl;
}