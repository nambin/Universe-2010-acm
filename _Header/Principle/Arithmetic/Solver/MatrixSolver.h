#ifndef _ARITHMETIC_MATRIX_SOLVER_H_
#define _ARITHMETIC_MATRIX_SOLVER_H_

#include "../../Principle.h"

namespace arithmetic {

class MatrixSolver {
public : 
	template<class Matrix, class VectorX, class VectorB, class Preconditioner, class Iteration, class VECOP>
	static void cg(const Matrix& A, VectorX& x, const VectorB& b, const Preconditioner& M, Iteration& iter, const VECOP& vecop);
};

template<class Matrix, class VectorX, class VectorB, class Preconditioner, class Iteration, class VECOP >
void MatrixSolver::cg(const Matrix& A, VectorX& x, const VectorB& b, const Preconditioner& M, Iteration& iter, const VECOP& vecop) {

	typedef VectorX TmpVec;
	typedef typename VectorX::value_type REAL_TYPE;
	typedef typename VectorX::index_type index_type;

	index_type size = x.size();
	REAL_TYPE rho(0), rho_new(0), alpha(0), beta(0);
	TmpVec s, r, z; s.resize(x); r.resize(x); z.resize(x);		
		
	mul(r, A, x);  vecop.dif(r, b, r);			// r = b - Ax	
	M.applyPreconditioner(s, r);				// s = M^-1 r
	rho = vecop.dot_directly(s, r);

	Timer timer;

	while(true) {

		mul(z, A, s);										// z = A * s;
		alpha = rho / vecop.dot_directly(s, z);
		
		vecop.icr_scaled_directly(x, s, +alpha);			// x += alpha * s;
		vecop.icr_scaled_directly(r, z, -alpha);			// r -= alpha * z;

		if(iter.finished(r)) break;

		M.applyPreconditioner(z, r);						// z = M^-1 r

		rho_new = vecop.dot_directly(z, r);
		beta = rho_new / rho;

		vecop.scaled_icr_directly(s, beta, z);			// s = z + beta * s	
		rho = rho_new;
	}

	timer.printTime("Time for CG Iteration");
	log_info << "Error for CG : " << iter.error << std::endl;
	iter.post_process();
}

};			// arithmetic

#endif