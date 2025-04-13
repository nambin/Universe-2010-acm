#ifndef _ARITHMETIC_GRID3_POISSON_SOLVER_H_
#define _ARITHMETIC_GRID3_POISSON_SOLVER_H_

#include "../../../Basis/ParallelProgramming/ParallelComputable.h"
#include "../../../Basis/Basis/Environments.h"
#include "../Solver/MatrixSolver.h"
#include "../Solver/Iteration.h"
#include "GridLibrary.h"
#include "Grid3_PoissonSolver_Functor.h"
#include "Grid3_PoissonMatrix.h"

namespace arithmetic {

template<typename T>
class Grid3_PoissonSolver : public ParallelComputable {
public : 
	typedef typename Grid3<T,T>::index_type index_type;

	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Grid3_PoissonSolver() {}
	Grid3_PoissonSolver(index_type nX, index_type nY, index_type nZ, T h, const Vector3<T>& offset) { resize(nX, nY, nZ, h, offset); }
	virtual ~Grid3_PoissonSolver() {}

	void setComputingOption(const ParallelComputable& rhs) { ParallelComputable::setComputingOption(rhs); poissonMatrix.setComputingOption(rhs); }

	index_type size() const { return x.size(); }
	
	void resize(index_type nX, index_type nY, index_type nZ, T h, const Vector3<T>& offset);
	template<typename T> void resize(const Grid3<T,T>& grid) { resize(grid.nX, grid.nY, grid.nZ, grid.dx, grid.offset); }

	virtual const bool importBinary(std::ifstream&) { return true; }
	virtual const bool exportBinary(std::ofstream&) const { return true; }

	//////////////////////////////////////////////////////////////////////////
	// Solver
	virtual const size_t solve(size_t maxItr, T tolerance) { basis::throwError("Not Implemented _ Grid3_PoissonSolver::solve()"); return 0; }

	const T normalize();
	virtual void scale(const T s);

public : 
	Grid3_PoissonMatrix<T>		poissonMatrix;
	Grid3<T,T>					x, b;
};

template<typename T>
const T Grid3_PoissonSolver<T>::normalize() {
	const T inf = poissonMatrix.poisson.infNorm();
	scale(1.0/inf);

	return 1.0 / inf;
}

template<typename T>
void Grid3_PoissonSolver<T>::scale(const T s) {
	poissonMatrix.poisson *= s;
	b *= s;
}

template<typename T>
void Grid3_PoissonSolver<T>::resize(index_type nX, index_type nY, index_type nZ, T h, const Vector3<T>& offset) {
	poissonMatrix.resize(nX,nY,nZ,h,offset);
	x.init(nX,nY,nZ,h,h,h,offset); b.init(nX,nY,nZ,h,h,h,offset);	
	x.zero(); b.zero();
}

};		// arithmetic

#endif