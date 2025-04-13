#pragma once

#include "Grid3_PoissonSolver.h"
#include "Grid3_Preconditioner_ICCG.h"
#include "Grid3_Preconditioner_BlockRedBlackICCG.h"	

namespace arithmetic {

template<typename T>
class Grid3_PoissonSolver_CG : public Grid3_PoissonSolver<T> {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Grid3_PoissonSolver_CG() { preconditioner = new Grid3_Preconditioner_Identity<T>(); }
	Grid3_PoissonSolver_CG(index_type nX, index_type nY, index_type nZ, T h, const Vector3<T>& offset) : Grid3_PoissonSolver<T>(nX, nY, nZ, h, offset) { preconditioner = new Grid3_Preconditioner_Identity<T>(); }
	~Grid3_PoissonSolver_CG() { if(preconditioner) delete preconditioner; }

	//////////////////////////////////////////////////////////////////////////
	// Make Preconditioner
	void setPreconditioner_ICCG();
	void setPreconditioner_BlockRedBlackICCG(index_type nbX, index_type nbY, index_type nbZ);

	//////////////////////////////////////////////////////////////////////////
	// Solver
	const size_t solve(size_t maxItr, T tolerance);

	static const std::string XmlElementName() { return std::string("Grid3_PoissonSolver_CG"); }

public : 
	Grid3_Preconditioner<T>*	preconditioner;
};

template<typename T>
void Grid3_PoissonSolver_CG<T>::setPreconditioner_ICCG() {	
	if(preconditioner) delete preconditioner;
	preconditioner = new Grid3_Preconditioner_ICCG<T>(poissonMatrix.nX(), poissonMatrix.nY(), poissonMatrix.nZ());
}

template<typename T>
void Grid3_PoissonSolver_CG<T>::setPreconditioner_BlockRedBlackICCG(index_type nbX, index_type nbY, index_type nbZ) {
	if(preconditioner) delete preconditioner;
	preconditioner = new Grid3_Preconditioner_BlockRedBlackICCG<T>(poissonMatrix.nX(), poissonMatrix.nY(), poissonMatrix.nZ(), nbX, nbY, nbZ);
}

template<typename T>
const size_t Grid3_PoissonSolver_CG<T>::solve(size_t maxItr, T tolerance) {

	VectorOperation<Grid3<T,T>> vecop(*this);
	Basic_Iteration<Grid3<T,T>> iter(maxItr, tolerance);
	iter.setComputingOption(*this);

	T rho(0), rho_new(0), alpha(0), beta(0);
	Grid3<T,T> s, r, z; s.resize(x); r.resize(x); z.resize(x);	
	
	preconditioner->buildPreconditioner(poissonMatrix);

	//	MatrixSolver::cg(poissonMatrix, x, b, *preconditioner, iter, VectorOperation<Grid3<T,T>>(*this));
	//	return iter.get_num_iter();

	mul(r, poissonMatrix, x);  
	vecop.dif(r, b, r);											// r = b - Ax	
	(*preconditioner).applyPreconditioner(s, r);				// s = M^-1 r

	rho = vecop.dot_directly(s, r);								// rho = s * r

	Timer timer;
	while(true) {

		T _dot = T(0.0);

#ifdef USE_OPENMP
		int my_end = int(s.nX-1);
		int i;
#pragma omp parallel for reduction(+:_dot) 
#else
		index_type my_end = s.nX-1;
		index_type i;
#endif
		for(i=1;i<my_end;++i) {
			for(index_type j=1;j<s.nY-1;++j) {
				for(index_type k=1;k<s.nZ-1;++k) {
					if( !poissonMatrix.isBoundary(i,j,k) )	{
						z(i,j,k)	= poissonMatrix.poisson(i,j,k,0)   * s(i,j,k)
									+ poissonMatrix.poisson(i,j,k-1,3) * s(i,j,k-1)
									+ poissonMatrix.poisson(i,j,k,3)   * s(i,j,k+1)
									+ poissonMatrix.poisson(i-1,j,k,1) * s(i-1,j,k)
									+ poissonMatrix.poisson(i,j,k,1)   * s(i+1,j,k)
									+ poissonMatrix.poisson(i,j-1,k,2) * s(i,j-1,k)
									+ poissonMatrix.poisson(i,j,k,2)   * s(i,j+1,k);		// z = A * s;
					}

					_dot += s(i,j,k) * z(i,j,k);			// _dot = s*z
				}}}

		alpha = rho / _dot;		

#ifdef USE_OPENMP
		my_end = int(s.nX-1);
#pragma omp parallel for 
#else
		my_end = s.nX-1;
#endif
		for(i=1;i<my_end;++i) {
			for(index_type j=1;j<s.nY-1;++j) {
				for(index_type k=1;k<s.nZ-1;++k) {
					x(i,j,k) += (+ alpha * s(i,j,k));		// x += alpha * s;
					r(i,j,k) += (- alpha * z(i,j,k));		// r -= alpha * z;
				}}}		

		if(iter.finished(r)) break;
		(*preconditioner).applyPreconditioner(z, r);		// z = M^-1 r

		rho_new = T(0.0);

#ifdef USE_OPENMP
		my_end = int(s.nX-1);
#pragma omp parallel for reduction(+:rho_new) 
#else
		my_end = s.nX-1;
#endif
		for(i=1;i<my_end;++i) {
			for(index_type j=1;j<s.nY-1;++j) {
				for(index_type k=1;k<s.nZ-1;++k) {
					rho_new += r(i,j,k) * z(i,j,k);			// rho_new = r * z
				}}}	
		beta = rho_new / rho;

#ifdef USE_OPENMP
		my_end = int(s.nX-1);
#pragma omp parallel for 
#else
		my_end = s.nX-1;
#endif
		for(i=1;i<my_end;++i) {
			for(index_type j=1;j<s.nY-1;++j) {
				for(index_type k=1;k<s.nZ-1;++k) {
					s(i,j,k) = beta * s(i,j,k) + z(i,j,k);	// s = z + beta * s	
				}}}

		rho = rho_new;
	}

	iter.post_process();
	log_info << "Error for Preconditioned CG : " << iter.error << std::endl;

	return iter.get_num_iter();
}

};		// arithmetic