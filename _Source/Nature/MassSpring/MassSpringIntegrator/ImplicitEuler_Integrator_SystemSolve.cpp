#include <stdafx.h>
#include <Nature/MassSpring/MassSpringIntegrator/ImplicitEuler_Integrator.h>
#include <Nature/MassSpring/Core/ConstrainedPCG.h>

using namespace std;
using namespace springspace;

void ImplicitEuler_Integrator::solveConstrainedPCG(const Constraints<real>& constraints) {
	ConstrainedPCG<real> pcgsolver(mat, constraints);
	pcgsolver.setParameters(1.0e-3, 100);

	real residual;
	size_t nIter;

	pcgsolver.pcgSolve(dv, b, residual, nIter);

	cout << "Residual : " << residual << endl;
	cout << "Iteration : " << nIter << endl;
}