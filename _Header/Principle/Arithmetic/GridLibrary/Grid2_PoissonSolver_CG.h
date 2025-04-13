#pragma once

#include "Grid2_PoissonSolver.h"
#include "Grid2_Preconditioner_ICCG.h"
#include "Grid2_Preconditioner_BlockRedBlackICCG.h"	

namespace arithmetic {

template<typename T>
class Grid2_PoissonSolver_CG : public Grid2_PoissonSolver<T> {	
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Grid2_PoissonSolver_CG() { preconditioner = new Grid2_Preconditioner_Identity<T>(); }
	Grid2_PoissonSolver_CG(index_type nX, index_type nY, T h, const Vector2<T>& offset) : Grid2_PoissonSolver<T>(nX, nY, h, offset) { preconditioner = new Grid2_Preconditioner_Identity<T>(); }
	~Grid2_PoissonSolver_CG() { if(preconditioner) delete preconditioner; }

	//////////////////////////////////////////////////////////////////////////
	// Make Preconditioner
	void setPreconditioner_ICCG();
	void setPreconditioner_BlockRedBlackICCG(index_type nbX, index_type nbY);

	//////////////////////////////////////////////////////////////////////////
	// Solver
	const size_t solve(size_t maxItr, T tolerance);	

	static const std::string XmlElementName() { return std::string("Grid2_PoissonSolver_CG"); }

public : 
	Grid2_Preconditioner<T>*	preconditioner;
};

template<typename T>
void Grid2_PoissonSolver_CG<T>::setPreconditioner_ICCG() {	
	if(preconditioner) delete preconditioner;
	preconditioner = new Grid2_Preconditioner_ICCG<T>(poissonMatrix.nX(), poissonMatrix.nY());
}

template<typename T>
void Grid2_PoissonSolver_CG<T>::setPreconditioner_BlockRedBlackICCG(index_type nbX, index_type nbY) {
	if(preconditioner) delete preconditioner;
	preconditioner = new Grid2_Preconditioner_BlockRedBlackICCG<T>(poissonMatrix.nX(), poissonMatrix.nY(), nbX, nbY);
}

// Need to Parallelize
template<typename T>
const size_t Grid2_PoissonSolver_CG<T>::solve(size_t maxItr, T tolerance) {
	Basic_Iteration<Grid2<T,T>> iteration(maxItr, tolerance);
	iteration.setComputingOption(*this);

	preconditioner->buildPreconditioner(poissonMatrix);
	MatrixSolver::cg(poissonMatrix, x, b, *preconditioner, iteration, VectorOperation<Grid2<T,T>>(*this));

	return iteration.get_num_iter();
}

};		// arithmetic