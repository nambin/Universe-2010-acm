#ifndef _ARITHMETIC_GRID2_POISSON_SOLVER_H_
#define _ARITHMETIC_GRID2_POISSON_SOLVER_H_

#include "../../../Basis/ParallelProgramming/ParallelComputable.h"
#include "../../../Basis/Basis/Environments.h"
#include "../Solver/MatrixSolver.h"
#include "../Solver/Iteration.h"
#include "GridLibrary.h"
#include "Grid2_PoissonSolver_Functor.h"
#include "Grid2_PoissonMatrix.h"

namespace arithmetic {

template<typename T>
class Grid2_PoissonSolver : public ParallelComputable {	
public : 
	typedef typename Grid2<T,T>::index_type index_type;

	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Grid2_PoissonSolver() {}
	Grid2_PoissonSolver(index_type nX, index_type nY, T h, const Vector2<T>& offset) { resize(nX, nY, h, offset); }
	virtual ~Grid2_PoissonSolver() {}

	void setComputingOption(const ParallelComputable& rhs) { ParallelComputable::setComputingOption(rhs); poissonMatrix.setComputingOption(rhs); }

	index_type size() const { return x.size(); }

	void resize(index_type nX, index_type nY, T h, const Vector2<T>& offset);
	template<typename R> void resize(const Grid2<R,T>& grid) { resize(grid.nX, grid.nY, grid.dx, grid.offset); }

	void set_initial_zero() { x.zero(); }

	virtual const T normalize();
	virtual void scale(const T s);

	virtual const bool importBinary(std::ifstream&);
	virtual const bool exportBinary(std::ofstream&) const;

	//////////////////////////////////////////////////////////////////////////
	// Solver
	virtual const size_t solve(size_t maxItr, T tolerance) { basis::throwError("Not Implemented _ Grid2_PoissonSolver::solve()"); return 0; }

public : 
	Grid2_PoissonMatrix<T>		poissonMatrix;
	Grid2<T,T>					x, b;
};

template<typename T>
const bool Grid2_PoissonSolver<T>::importBinary(std::ifstream& ifs) {
	poissonMatrix.importBinary(ifs);
	b.importBinary(ifs);
	x.resize(b);

	return true;
}

template<typename T>
const bool Grid2_PoissonSolver<T>::exportBinary(std::ofstream& ofs) const {
	poissonMatrix.exportBinary(ofs);
	b.exportBinary(ofs);

	return true;
}

template<typename T>
const T Grid2_PoissonSolver<T>::normalize() {
	const T inf = poissonMatrix.poisson.infNorm();
	scale(1.0/inf);

	return 1.0 / inf;
}

template<typename T>
void Grid2_PoissonSolver<T>::scale(const T s) {
	poissonMatrix.poisson *= s;
	b *= s;
}

template<typename T>
void Grid2_PoissonSolver<T>::resize(index_type nX, index_type nY, T h, const Vector2<T>& offset) {
	poissonMatrix.resize(nX, nY, h, offset);
	x.init(nX,nY,h,h,offset); b.init(nX,nY,h,h,offset);	
	x.zero(); b.zero();
}

};		// arithmetic

#endif