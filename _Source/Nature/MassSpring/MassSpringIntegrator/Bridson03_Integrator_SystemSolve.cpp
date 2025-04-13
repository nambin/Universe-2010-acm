#include <stdafx.h>
#include <Nature/MassSpring/MassSpringIntegrator/Bridson03_Integrator.h>

using namespace std;
using namespace springspace;

bool Bridson03_Integrator::pcgSolve(Vec3Array<real>& result, const Vec3Array<real>& rhs, const TriDynamicMesh<real>& mesh, real dt, real &residual_out, size_t &iterations_out) {
	const Constraints<real>& constraints = mesh.constraint;

	size_t n = b.size();

	for(size_t i=0;i<n;++i)	
		result[i] = constraints.getVelocity(i);						// dv = constraint vel

	Vector3<real> tmp;
	residual_out = 0;
	for(size_t i=0;i<n;++i)	{
		constraints.filterVelocity(i, tmp = rhs[i]);
		residual_out += (arithmetic::sqrmag(tmp) / invDiag[i]);
	}	
	//cout << "Initial Residual : " << residual_out << endl;

	mul(z, mesh, dt, result);
	for(size_t i=0;i<n;++i)	
		r[i] = rhs[i] - z[i];
	constraints.filterVelocities(r);								// r = filter(b - mat * dv)

	applyDiagonalPreconditioner(s, r);
	constraints.filterVelocities(s);								// s = filter(P^-1 r)

	real rho = dot(s, r);
	if(rho == 0) {
		iterations_out = 0;
		return true;
	}

	real tol = tolerance * tolerance * residual_out;

	size_t iteration;
	real alpha, beta, rho_new;

	for(iteration=0; iteration < maxIter; ++iteration) {
		residual_out = rho;
		if(residual_out <= tol) {
			iterations_out = iteration+1;
			return true; 
		}

		mul(z, mesh, dt, s); constraints.filterVelocities(z);	// z = filter(mat * s)
		alpha = rho / dot(s, z);								// alpha = rho / dot(s, z)
		increment(result, alpha, s);							// result += s * alpha;
		increment(r, -alpha, z);								// r -= z * alpha;

		applyDiagonalPreconditioner(z, r);						// z = P^-1 * r
		rho_new = dot(z, r);
		beta = rho_new / rho;
		s *= beta; s += z; constraints.filterVelocities(s);		// s = filter(z+beta*s)
		rho = rho_new;
	}

	iterations_out = iteration;
	return false;
}

void Bridson03_Integrator::makeDiagonalPreconditioner(const TriDynamicMesh<real>& mesh, real dt) {
	real damp = mesh.stretchFunc->dampingCoeff();
	real tmp = damp * dt; 

	for(size_t i=0;i<invDiag.size();++i) {
		invDiag[i] = mesh.mass[i] + tmp * mesh.s2v[i].size();
		invDiag[i] = 1.0f / invDiag[i];
	}
}

void Bridson03_Integrator::applyDiagonalPreconditioner(Vec3Array<real>& ret, const Vec3Array<real>& b) const {
	for(size_t i=0;i<ret.size();++i)
		ret[i] = b[i] * invDiag[i];
}

void Bridson03_Integrator::mul(Vec3Array<real> &ret, const TriDynamicMesh<real>& mesh, real dt, const Vec3Array<real> &v) {
	real damp = mesh.stretchFunc->dampingCoeff();
	real tmp = damp * dt; 

	for(size_t i=0;i<ret.size();++i) {
		ret[i] = v[i] * mesh.mass[i];

		for(size_t k=0;k<mesh.v2v[i].size();++k) {
			ret[i] += (tmp * (v[i] - v[mesh.v2v[i][k]]));
		}
	}
}