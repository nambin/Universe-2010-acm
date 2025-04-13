#ifndef _MASS_SPRING_CONSTRAINED_PCG_H_
#define _MASS_SPRING_CONSTRAINED_PCG_H_

#include "../Core/MassSpringSystemMatrix.h"

namespace springspace {

template<typename T>
class DiagonalPreconditioner {
public : 
	explicit DiagonalPreconditioner(size_t _n = 0);

	void clear() { resize(0); }
	void resize(size_t);
	
	void makePreconditioner(const MassSpringSystemMatrix<T>&);
	void applyPreconditioner(Vec3Array<T>& ret, const Vec3Array<T>& b) const;

public : 
	size_t			n;
	Vec3Array<T>	invDiag;
};

template<typename T>
class ConstrainedPCG {
public : 
	ConstrainedPCG(const MassSpringSystemMatrix<T>&, const Constraints<T>&);
	void allocateMemory();
	void setParameters(T tol, size_t maxIter);

	bool pcgSolve(Vec3Array<T>& result, const Vec3Array<T>& rhs, T &residual_out, size_t &iterations_out);

protected:
	const MassSpringSystemMatrix<T>&	mat;				// used within loop
	const Constraints<T>&				constraints;

	DiagonalPreconditioner<T>			diagPreconditioner;
	Vec3Array<T>						z, s, r;			// temporary vectors for PCG	

	// parameters
	T		tolerance;
	size_t	maxIter;
};

template<typename T>
DiagonalPreconditioner<T>::DiagonalPreconditioner(size_t _n) { resize(_n); }

template<typename T>
void DiagonalPreconditioner<T>::resize(size_t _n) {
	n = _n;	invDiag.resize(n);
}

template<typename T>
void DiagonalPreconditioner<T>::makePreconditioner(const MassSpringSystemMatrix<T>& mat) {

	resize(mat.n);
	for(size_t i=0;i<n;++i) {
		invDiag[i].zero();
		for(Matrix3<T>::index_type k=0;k<3;++k) {
			invDiag[i](k) = 1.0f / mat.diag[i](k,k);
		}
	}

	//T _min = 10.0e+5;
	//for(size_t i=0;i<n;++i) 
	//	for(size_t k=0;k<3;++k) 
	//		_min = (_min < mat.diag[i](k,k)) ? _min : mat.diag[i](k,k);
	//cout << "Minimum of Diagonal : " << _min << endl;
}

template<typename T>
void DiagonalPreconditioner<T>::applyPreconditioner(Vec3Array<T>& ret, const Vec3Array<T>& b) const {
	for(size_t i=0;i<n;++i) {
		ret[i][0] = invDiag[i][0] * b[i][0];
		ret[i][1] = invDiag[i][1] * b[i][1];
		ret[i][2] = invDiag[i][2] * b[i][2];
	}
}

template<typename T>
ConstrainedPCG<T>::ConstrainedPCG(const MassSpringSystemMatrix<T>& _mat, const Constraints<T>& _ctrs)
: mat(_mat), constraints(_ctrs)
{
	diagPreconditioner.makePreconditioner(mat);
	allocateMemory();
	setParameters(1.0e-5, 300);
}

template<typename T>
void ConstrainedPCG<T>::allocateMemory() {
	z.resize(mat.n); s.resize(mat.n); r.resize(mat.n);
}

template<typename T>
void ConstrainedPCG<T>::setParameters(T _tol, size_t _maxIter) {
	tolerance = _tol; maxIter = _maxIter;
}

template<typename T>
bool ConstrainedPCG<T>::pcgSolve(Vec3Array<T>& result, const Vec3Array<T>& rhs, T &residual_out, size_t &iterations_out) {
	size_t n = mat.n;

	for(size_t i=0;i<n;++i)	
		result[i] = constraints.getVelocity(i);						// dv = constraint vel

	Vector3<T> tmp;
	residual_out = 0;
	for(size_t i=0;i<n;++i)	{
		constraints.filterVelocity(i, tmp = rhs[i]);
		for(Vector3<T>::index_type k=0;k<3;++k)
			residual_out += tmp[k] * tmp[k] * mat.diag[i](k,k);		// residual_out = filter(rhs)^T * P * filter(rhs)
	}	
//	cout << "Initial Residual : " << residual_out << endl;

	mul(z, mat, result);
	for(size_t i=0;i<n;++i)	
		r[i] = rhs[i] - z[i];
	constraints.filterVelocities(r);								// r = filter(b - mat * dv)

	diagPreconditioner.applyPreconditioner(s, r);
	constraints.filterVelocities(s);								// s = filter(P^-1 r)
	
	T rho = dot(s, r);
	if(rho == 0) {
		iterations_out = 0;
		return true;
	}

	T tol = tolerance * tolerance * residual_out;

	size_t iteration;
	T alpha, beta, rho_new;

	for(iteration=0; iteration < maxIter; ++iteration) {
		residual_out = rho;
		if(residual_out <= tol) {
			iterations_out = iteration+1;
			return true; 
		}

		mul(z, mat, s);	constraints.filterVelocities(z);		// z = filter(mat * s)
		alpha = rho / dot(s, z);								// alpha = rho / dot(s, z)
		increment(result, alpha, s);							// result += s * alpha;
		increment(r, -alpha, z);								// r -= z * alpha;

		diagPreconditioner.applyPreconditioner(z, r);			// z = P^-1 * r
		rho_new = dot(z, r);
		beta = rho_new / rho;
		s *= beta; s += z; constraints.filterVelocities(s);		// s = filter(z+beta*s)
		rho = rho_new;
	}

	iterations_out = iteration;
	return false;
}
	
};			// springspace

#endif