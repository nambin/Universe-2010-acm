#include <stdafx.h>
#include <Nature/Fluid/CFDsolver/CFDProjection.h>

using namespace std;
using namespace fluidspace;

void CFDProjection_SinglePhase::buildPoissonMatrix(Grid2_PoissonSolver<double>& solver, const CFDsolver2D& cfd) {

	Grid2d u_part, v_part;
	buildPoissonMatrix(u_part, v_part, cfd);

	if(TypeChecker<Grid2_PoissonSolver_CG<double>>()(&solver)) 
		CFDProjection::buildPoissonMatrix(solver.poissonMatrix.poisson, solver.poissonMatrix.isBoundary, u_part, v_part);

	else if(TypeChecker<Grid2_PoissonSolver_GMG<double>>()(&solver)) 
		CFDProjection::buildPoissonMatrix(static_cast<Grid2_PoissonSolver_GMG<double>&>(solver), u_part, v_part);	
}

void CFDProjection_SinglePhase::buildPoissonMatrix(Grid3_PoissonSolver<double>& solver, const CFDsolver3D& cfd) {

	Grid3d u_part, v_part, w_part;
	buildPoissonMatrix(u_part, v_part, w_part, cfd);

	if(TypeChecker<Grid3_PoissonSolver_CG<double>>()(&solver)) 
		CFDProjection::buildPoissonMatrix(solver.poissonMatrix.poisson, solver.poissonMatrix.isBoundary, u_part, v_part, w_part);

	else if(TypeChecker<Grid3_PoissonSolver_GMG<double>>()(&solver)) 
		CFDProjection::buildPoissonMatrix(static_cast<Grid3_PoissonSolver_GMG<double>&>(solver), u_part, v_part, w_part);
}

void CFDProjection_SinglePhase::buildPoissonMatrix(Grid2d& u_part, Grid2d& v_part, const CFDsolver2D& cfd) {
	if(u_part.is_empty())	u_part.resize(cfd.u);
	if(v_part.is_empty())	v_part.resize(cfd.v);

	const Grid2b&	isKinematicBoundary = cfd.isKinematicBoundary;
	const double hInv_sq = cfd.hInv * cfd.hInv;
	const double one = 1.0;

	u_part.zero();
	v_part.zero();

	for( index_type i=0; i < isKinematicBoundary.nX ; ++i ) {
		for( index_type j=0; j < isKinematicBoundary.nY; ++j ) {
			if( isKinematicBoundary(i,j) == false) {
				//if( !isKinematicBoundary(i+1,j) )	u_part(i,j) = hInv_sq;
				//if( !isKinematicBoundary(i,j+1) )	v_part(i,j) = hInv_sq;

				if( !isKinematicBoundary(i+1,j) )	u_part(i,j) = one;
				if( !isKinematicBoundary(i,j+1) )	v_part(i,j) = one;
			}
		}
	}	
}

void CFDProjection_SinglePhase::buildPoissonMatrix(Grid3d& u_part, Grid3d& v_part, Grid3d& w_part, const CFDsolver3D& cfd) {
	if(u_part.is_empty())	u_part.resize(cfd.u);
	if(v_part.is_empty())	v_part.resize(cfd.v);
	if(w_part.is_empty())	w_part.resize(cfd.w);

	const Grid3b&	isKinematicBoundary = cfd.isKinematicBoundary;
	const double hInv_sq = cfd.hInv * cfd.hInv;
	const double one = 1.0;

	u_part.zero();
	v_part.zero();
	w_part.zero();

#ifdef USE_OPENMP
	int my_end = int(isKinematicBoundary.nX);
	int i;

#pragma omp parallel for
#else
	index_type my_end = isKinematicBoundary.nX;
	index_type i;
#endif
	for( i = 0; i < my_end; ++i ) {	
		for( index_type j=0; j < isKinematicBoundary.nY; ++j ) {
			for( index_type k=0; k < isKinematicBoundary.nZ; ++k ) {
				if( isKinematicBoundary(i,j,k) == false) {
					//if( !isKinematicBoundary(i+1,j,k) )	u_part(i,j,k) = hInv_sq;
					//if( !isKinematicBoundary(i,j+1,k) )	v_part(i,j,k) = hInv_sq;
					//if( !isKinematicBoundary(i,j,k+1) )	w_part(i,j,k) = hInv_sq;

					if( !isKinematicBoundary(i+1,j,k) )	u_part(i,j,k) = one;
					if( !isKinematicBoundary(i,j+1,k) )	v_part(i,j,k) = one;
					if( !isKinematicBoundary(i,j,k+1) )	w_part(i,j,k) = one;
				}
			}
		}
	}	
}