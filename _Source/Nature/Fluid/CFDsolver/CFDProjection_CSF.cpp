#include <stdafx.h>
#include <Nature/Fluid/CFDsolver/CFDProjection.h>
#include <Principle/Geometry/MacGrid/MacGridLibrary.h>

using namespace std;
using namespace fluidspace;

void CFDProjection_CSF::buildPoissonMatrix(Grid2_PoissonSolver<double>& solver, const Grid2r& phi, const CFDsolver2D& cfd) {
	Grid2d u_part, v_part;
	buildPoissonMatrix(u_part, v_part, phi, cfd);

	if(TypeChecker<Grid2_PoissonSolver_CG<double>>()(&solver)) 
		CFDProjection::buildPoissonMatrix(solver.poissonMatrix.poisson, solver.poissonMatrix.isBoundary, u_part, v_part);

	else if(TypeChecker<Grid2_PoissonSolver_GMG<double>>()(&solver)) 
		CFDProjection::buildPoissonMatrix(static_cast<Grid2_PoissonSolver_GMG<double>&>(solver), u_part, v_part);	
}

void CFDProjection_CSF::buildPoissonMatrix(Grid3_PoissonSolver<double>& solver, const Grid3r& phi, const CFDsolver3D& cfd) {
	Grid3d u_part, v_part, w_part;
	buildPoissonMatrix(u_part, v_part, w_part, phi, cfd);

	if(TypeChecker<Grid3_PoissonSolver_GMG<double>>()(&solver)) 
		CFDProjection::buildPoissonMatrix(static_cast<Grid3_PoissonSolver_GMG<double>&>(solver), u_part, v_part, w_part);	

	else if(TypeChecker<Grid3_PoissonSolver_CG<double>>()(&solver)) 
		CFDProjection::buildPoissonMatrix(solver.poissonMatrix.poisson, solver.poissonMatrix.isBoundary, u_part, v_part, w_part);
}

void CFDProjection_CSF::buildPoissonMatrix(Grid2d& u_part, Grid2d& v_part, const Grid2r& phi, const CFDsolver2D& cfd) {

	const Grid2b& isKinematicBoundary = cfd.isKinematicBoundary;

	const real hInv = cfd.hInv;
	const real hInv_sq = cfd.hInv * cfd.hInv;	

	if(u_part.is_empty())	u_part.resize(cfd.u);
	if(v_part.is_empty())	v_part.resize(cfd.v);

	u_part.zero();
	v_part.zero();

	for( index_type i=1; i < isKinematicBoundary.nX-1 ; ++i ) {
		real _density, _phi;

		for( index_type j=1; j < isKinematicBoundary.nY-1 ; ++j ) {
			if( isKinematicBoundary(i,j) == false) {	
				if( !isKinematicBoundary(i+1,j) ) {
					_phi = 0.5f * (phi(i,j) + phi(i+1,j));
					_density = MacGridLibrary<real>::heaviside(_phi * hInv, rho_pos, rho_neg);		

					//u_part(i,j) = (1.0 / _density) * hInv_sq;
					u_part(i,j) = (1.0 / _density);
				}
				if( !isKinematicBoundary(i,j+1) ) {
					_phi = 0.5f * (phi(i,j) + phi(i,j+1));
					_density = MacGridLibrary<real>::heaviside(_phi * hInv, rho_pos, rho_neg);		

					//v_part(i,j) = (1.0 / _density) * hInv_sq;
					v_part(i,j) = (1.0 / _density);
				}
			}
		}
	}	
}

void CFDProjection_CSF::buildPoissonMatrix(Grid3d& u_part, Grid3d& v_part, Grid3d& w_part, const Grid3r& phi, const CFDsolver3D& cfd) {
	const Grid3b& isKinematicBoundary = cfd.isKinematicBoundary;

	const real h_half = cfd.h * 0.5;
	const real hInv = cfd.hInv;
	const real hInv_sq = cfd.hInv * cfd.hInv;	

	if(u_part.is_empty())	u_part.resize(cfd.u);
	if(v_part.is_empty())	v_part.resize(cfd.v);
	if(w_part.is_empty())	w_part.resize(cfd.w);

	u_part.zero();
	v_part.zero();
	w_part.zero();

	for( index_type i=1; i < isKinematicBoundary.nX-1 ; ++i ) {
		real _density, _phi;

		for( index_type j=1; j < isKinematicBoundary.nY-1 ; ++j ) {
			for( index_type k=1; k < isKinematicBoundary.nZ-1 ; ++k ) {

				if( isKinematicBoundary(i,j,k) == false) {	
					if( !isKinematicBoundary(i+1,j,k) ) {
						_phi = 0.5f * (phi(i,j,k) + phi(i+1,j,k));
						_density = MacGridLibrary<real>::heaviside(_phi * hInv, rho_pos, rho_neg);		

						//u_part(i,j,k) = (1.0 / _density) * hInv_sq;
						u_part(i,j,k) = (1.0 / _density) ;
					}
					if( !isKinematicBoundary(i,j+1,k) ) {
						_phi = 0.5f * (phi(i,j,k) + phi(i,j+1,k));
						_density = MacGridLibrary<real>::heaviside(_phi * hInv, rho_pos, rho_neg);		

						//v_part(i,j,k) = (1.0 / _density) * hInv_sq;
						v_part(i,j,k) = (1.0 / _density) ;
					}
					if( !isKinematicBoundary(i,j,k+1) ) {
						_phi = 0.5f * (phi(i,j,k) + phi(i,j,k+1));
						_density = MacGridLibrary<real>::heaviside(_phi * hInv, rho_pos, rho_neg);		

						//w_part(i,j,k) = (1.0 / _density) * hInv_sq;
						w_part(i,j,k) = (1.0 / _density) ;
					}
				}
			}
		}
	}	
}

void CFDProjection_CSF::addPressureGradient(CFDsolver2D& cfd, const Grid2d& pressure, const Grid2r& phi) {
	const Grid2b&		isKinematicBoundary = cfd.isKinematicBoundary;
	Grid2r&				u = cfd.u;
	Grid2r&				v = cfd.v;	

	const real hInv = cfd.hInv;

	for ( int i = 1; i < u.nX-1; ++i ) {
		real _density, _phi;

		for ( int j = 1; j < u.nY-1; ++j ) {		
			if( !isKinematicBoundary[cfd.iU2S(i,j,1)] && !isKinematicBoundary[cfd.iU2S(i,j,0)] )	{
				_phi = 0.5f * (phi(cfd.iU2S(i,j,0)) + phi(cfd.iU2S(i,j,1)));
				_density = MacGridLibrary<real>::heaviside(_phi*hInv, rho_pos, rho_neg);

				//u(i,j) += ( pressure[cfd.iU2S(i,j,1)] - pressure[cfd.iU2S(i,j,0)] ) / _density * hInv;
				u(i,j) += ( pressure[cfd.iU2S(i,j,1)] - pressure[cfd.iU2S(i,j,0)] ) / _density;
			}
		}
	}

	for ( int i = 1; i < v.nX-1; ++i ) {
		real _density, _phi;

		for ( int j = 1; j < v.nY-1; ++j ) {
			if( !isKinematicBoundary[cfd.iV2S(i,j,1)] && !isKinematicBoundary[cfd.iV2S(i,j,0)] )	{
				_phi = 0.5f * (phi(cfd.iV2S(i,j,0)) + phi(cfd.iV2S(i,j,1)));
				_density = MacGridLibrary<real>::heaviside(_phi*hInv, rho_pos, rho_neg);

				//v(i,j) += ( pressure[cfd.iV2S(i,j,1)] - pressure[cfd.iV2S(i,j,0)] ) / _density * hInv;
				v(i,j) += ( pressure[cfd.iV2S(i,j,1)] - pressure[cfd.iV2S(i,j,0)] ) / _density;
			}
		}
	}
}

void CFDProjection_CSF::addPressureGradient(CFDsolver3D& cfd, const Grid3d& pressure, const Grid3r& phi) {
	const Grid3b&		isKinematicBoundary = cfd.isKinematicBoundary;
	Grid3r&				u = cfd.u;
	Grid3r&				v = cfd.v;	
	Grid3r&				w = cfd.w;	

	const real h_half = cfd.h * 0.5;
	const real hInv = cfd.hInv;

	for ( int i = 1; i < u.nX-1; ++i ) {
		real _density, _phi;

		for ( int j = 1; j < u.nY-1; ++j ) {	
			for ( int k = 1; k < u.nZ-1; ++k ) {	
				if( !isKinematicBoundary[cfd.iU2S(i,j,k,1)] && !isKinematicBoundary[cfd.iU2S(i,j,k,0)] )	{
					_phi = 0.5f * (phi(cfd.iU2S(i,j,k,0)) + phi(cfd.iU2S(i,j,k,1)));
					_density = MacGridLibrary<real>::heaviside(_phi*hInv, rho_pos, rho_neg);

					//u(i,j,k) += ( pressure[cfd.iU2S(i,j,k,1)] - pressure[cfd.iU2S(i,j,k,0)] ) / _density * hInv;
					u(i,j,k) += ( pressure[cfd.iU2S(i,j,k,1)] - pressure[cfd.iU2S(i,j,k,0)] ) / _density;
				}
			}			
		}
	}

	for ( int i = 1; i < v.nX-1; ++i ) {
		real _density, _phi;

		for ( int j = 1; j < v.nY-1; ++j ) {
			for ( int k = 1; k < v.nZ-1; ++k ) {
				if( !isKinematicBoundary[cfd.iV2S(i,j,k,1)] && !isKinematicBoundary[cfd.iV2S(i,j,k,0)] )	{
					_phi = 0.5f * (phi(cfd.iV2S(i,j,k,0)) + phi(cfd.iV2S(i,j,k,1)));
					_density = MacGridLibrary<real>::heaviside(_phi*hInv, rho_pos, rho_neg);

					//v(i,j,k) += ( pressure[cfd.iV2S(i,j,k,1)] - pressure[cfd.iV2S(i,j,k,0)] ) / _density * hInv;
					v(i,j,k) += ( pressure[cfd.iV2S(i,j,k,1)] - pressure[cfd.iV2S(i,j,k,0)] ) / _density;
				}
			}
		}
	}

	for ( int i = 1; i < w.nX-1; ++i ) {
		real _density, _phi;

		for ( int j = 1; j < w.nY-1; ++j ) {
			for ( int k = 1; k < w.nZ-1; ++k ) {
				if( !isKinematicBoundary[cfd.iW2S(i,j,k,1)] && !isKinematicBoundary[cfd.iW2S(i,j,k,0)] )	{
					_phi = 0.5f * (phi(cfd.iW2S(i,j,k,0)) + phi(cfd.iW2S(i,j,k,1)));
					_density = MacGridLibrary<real>::heaviside(_phi*hInv, rho_pos, rho_neg);

					//w(i,j,k) += ( pressure[cfd.iW2S(i,j,k,1)] - pressure[cfd.iW2S(i,j,k,0)] ) / _density * hInv;
					w(i,j,k) += ( pressure[cfd.iW2S(i,j,k,1)] - pressure[cfd.iW2S(i,j,k,0)] ) / _density;
				}
			}
		}
	}
}