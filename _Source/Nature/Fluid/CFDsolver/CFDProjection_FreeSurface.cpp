#include <stdafx.h>
#include <Nature/Fluid/CFDsolver/CFDProjection.h>

using namespace std;
using namespace fluidspace;

void CFDProjection_FreeSurface::set_initial_guess(Grid2_PoissonSolver<double>& sol) { 
	//sol.x.zero(); return;

	const double zero(0.0);

	for(index_type i=1;i<sol.poissonMatrix.isBoundary.nX-1;++i) {
		for(index_type j=1;j<sol.poissonMatrix.isBoundary.nY-1;++j) {
			if(sol.poissonMatrix.isBoundary(i,j) == true)
				sol.x(i,j) = zero;
		}
	}
}

void CFDProjection_FreeSurface::set_initial_guess(Grid3_PoissonSolver<double>& sol) {
	//sol.x.zero(); return;

	const double zero(0.0);

	for(index_type i=1;i<sol.poissonMatrix.isBoundary.nX-1;++i) {
		for(index_type j=1;j<sol.poissonMatrix.isBoundary.nY-1;++j) {
			for(index_type k=1;k<sol.poissonMatrix.isBoundary.nZ-1;++k) {
				if(sol.poissonMatrix.isBoundary(i,j,k) == true)
					sol.x(i,j,k) = zero;
			}			
		}
	}
}

void CFDProjection_FreeSurface::buildDivergence(Grid2d& divergence, const Grid2r& phi, const CFDsolver2D& cfd) {
		
	const double hInv = divergence.dxInv;
	const real zero = real(0.0);

	divergence.zero();

	for ( index_type i = 1; i < divergence.nX-1; ++i ) {

		double rurx, rvry;
		Vector2<real> pos;

		for ( index_type j = 1; j < divergence.nY-1; ++j ) {

			pos = Vector2<real>(divergence.pos(i,j));

			if(phi(pos) < zero && !cfd.isKinematicBoundary(i,j)) {
				rurx = ( cfd.u(cfd.iS2U(i,j,1)) - cfd.u(cfd.iS2U(i,j,0)) );
				rvry = ( cfd.v(cfd.iS2V(i,j,1)) - cfd.v(cfd.iS2V(i,j,0)) );

				divergence(i,j) = (rurx + rvry);
				//divergence(i,j) = (rurx + rvry) * hInv;
			}
		}
	}
}

void CFDProjection_FreeSurface::buildPoissonMatrix(Grid2d& u_part, Grid2d& v_part, Grid2d& s_part, const Grid2r& phi, const CFDsolver2D& cfd) {

	const real zero = real(0.0);

	if(u_part.is_empty())	u_part.resize(cfd.u);
	if(v_part.is_empty())	v_part.resize(cfd.v);
	if(s_part.is_empty())	s_part.resize(cfd.isKinematicBoundary);

	const Grid2b&	isKinematicBoundary = cfd.isKinematicBoundary;
	const double	hInv_sq = cfd.hInv * cfd.hInv;
	const double	one = 1.0;

	u_part.zero();
	v_part.zero();
	s_part.zero();

	for ( index_type i = 1; i < isKinematicBoundary.nX-1; ++i ) {
		for ( index_type j = 1; j < isKinematicBoundary.nY-1; ++j ) {
			if( !isKinematicBoundary(i,j) && phi(i,j) < zero )	{				
				if( !isKinematicBoundary(i-1,j) && phi(i-1,j) >= zero )	s_part(i,j) += one;

				if( !isKinematicBoundary(i+1,j) && phi(i+1,j) >= zero )	s_part(i,j) += one;
				else if( !isKinematicBoundary(i+1,j) )					u_part(i,j) += one;

				if( !isKinematicBoundary(i,j-1) && phi(i,j-1) >= zero )	s_part(i,j) += one;		

				if( !isKinematicBoundary(i,j+1) && phi(i,j+1) >= zero )	s_part(i,j) += one;
				else if( !isKinematicBoundary(i,j+1) )					v_part(i,j) += one;				
			}
		}
	}	
}

void CFDProjection_FreeSurface::buildPoissonMatrix(Grid2_PoissonSolver<double>& solver, const Grid2r& phi, const CFDsolver2D& cfd) {
	Grid2d u_part, v_part, s_part;
	buildPoissonMatrix(u_part, v_part, s_part, phi, cfd);

	if(TypeChecker<Grid2_PoissonSolver_CG<double>>()(&solver)) 
		CFDProjection::buildPoissonMatrix(solver.poissonMatrix.poisson, solver.poissonMatrix.isBoundary, u_part, v_part, s_part);

	else if(TypeChecker<Grid2_PoissonSolver_GMG<double>>()(&solver)) 
		CFDProjection::buildPoissonMatrix(static_cast<Grid2_PoissonSolver_GMG<double>&>(solver), u_part, v_part, s_part);
}

void CFDProjection_FreeSurface::buildDivergence(Grid3d& divergence, const Grid3r& phi, const CFDsolver3D& cfd) {

	const real zero = real(0.0);

	divergence.zero();

#ifdef USE_OPENMP
	int my_end = int(divergence.nX-1);
	int i;

#pragma omp parallel for
#else
	index_type my_end = divergence.nX-1;
	index_type i;
#endif
	for( i = 1; i < my_end; ++i ) {

		double rurx, rvry, rwrz;
		for ( index_type j = 1; j < divergence.nY-1; ++j ) {
			for( index_type k = 1; k < divergence.nZ-1; ++k)	{
				if(phi(i,j,k) < zero && !cfd.isKinematicBoundary(i,j,k))	{
					rurx = ( cfd.u[cfd.iS2U(i,j,k,1)] - cfd.u[cfd.iS2U(i,j,k,0)] );
					rvry = ( cfd.v[cfd.iS2V(i,j,k,1)] - cfd.v[cfd.iS2V(i,j,k,0)] );
					rwrz = ( cfd.w[cfd.iS2W(i,j,k,1)] - cfd.w[cfd.iS2W(i,j,k,0)] );

					//divergence(i,j,k) = (rurx + rvry + rwrz) * cfd.hInv;	
					divergence(i,j,k) = (rurx + rvry + rwrz) ;	
				}
			}			
		}
	}
}

void CFDProjection_FreeSurface::buildPoissonMatrix(Grid3d& u_part, Grid3d& v_part, Grid3d& w_part, Grid3d& s_part, const Grid3r& phi, const CFDsolver3D& cfd) {

	if(u_part.is_empty())	u_part.resize(cfd.u);
	if(v_part.is_empty())	v_part.resize(cfd.v);
	if(w_part.is_empty())	w_part.resize(cfd.w);
	if(s_part.is_empty())	s_part.resize(cfd.isKinematicBoundary);

	const Grid3b&	isKinematicBoundary = cfd.isKinematicBoundary;
	const double	hInv_sq = cfd.hInv * cfd.hInv;
	const double	one = 1.0;

	const real		zero = real(0.0);	

	u_part.zero();
	v_part.zero();
	w_part.zero();
	s_part.zero();

#ifdef USE_OPENMP
	int my_end = int(isKinematicBoundary.nX-1);
	int i;

#pragma omp parallel for
#else
	index_type my_end = isKinematicBoundary.nX-1;
	index_type i;
#endif
	for( i = 1; i < my_end; ++i ) {
		for ( index_type j = 1; j < isKinematicBoundary.nY-1; ++j ) {
			for ( index_type k = 1; k < isKinematicBoundary.nZ-1; ++k ) {

				if( !isKinematicBoundary(i,j,k) && phi(i,j,k) < zero )	{				
					if( !isKinematicBoundary(i-1,j,k) && phi(i-1,j,k) >= zero )	s_part(i,j,k) += one;

					if( !isKinematicBoundary(i+1,j,k) && phi(i+1,j,k) >= zero )	s_part(i,j,k) += one;
					else if( !isKinematicBoundary(i+1,j,k) )					u_part(i,j,k) += one;

					if( !isKinematicBoundary(i,j-1,k) && phi(i,j-1,k) >= zero )	s_part(i,j,k) += one;		

					if( !isKinematicBoundary(i,j+1,k) && phi(i,j+1,k) >= zero )	s_part(i,j,k) += one;
					else if( !isKinematicBoundary(i,j+1,k) )					v_part(i,j,k) += one;

					if( !isKinematicBoundary(i,j,k-1) && phi(i,j,k-1) >= zero )	s_part(i,j,k) += one;

					if( !isKinematicBoundary(i,j,k+1) && phi(i,j,k+1) >= zero )	s_part(i,j,k) += one;
					else if( !isKinematicBoundary(i,j,k+1) )					w_part(i,j,k) += one;
				}

			}
		}
	}	
}

void CFDProjection_FreeSurface::buildPoissonMatrix(Grid3_PoissonSolver<double>& solver, const Grid3r& phi, const CFDsolver3D& cfd) {
	Grid3d u_part, v_part, w_part, s_part;
	buildPoissonMatrix(u_part, v_part, w_part, s_part, phi, cfd);
	
	if(TypeChecker<Grid3_PoissonSolver_CG<double>>()(&solver)) 
		CFDProjection::buildPoissonMatrix(solver.poissonMatrix.poisson, solver.poissonMatrix.isBoundary, u_part, v_part, w_part, s_part);

	else if(TypeChecker<Grid3_PoissonSolver_GMG<double>>()(&solver)) 
		CFDProjection::buildPoissonMatrix(static_cast<Grid3_PoissonSolver_GMG<double>&>(solver), u_part, v_part, w_part, s_part);
}

//void CFDProjection_FreeSurface::buildPoissonMatrix(Grid2_PoissonMatrix<double>& matrix, const Grid2r& phi, const CFDsolver2D& cfd) {
//
//	const real				zero = real(0.0);
//
//	Grid2x3d&				poisson = matrix.poisson;
//	Grid2<bool,double>&		isBoundary = matrix.isBoundary;
//
//	const AABB2D<real> aabb = phi.getAABB();
//
//	isBoundary.set_boundary(true);
//	for(index_type i=1;i<isBoundary.nX-1;++i) {
//		Vector2<real> pos;
//		for(index_type j=1;j<isBoundary.nY-1;++j) {
//			pos = Vector2<real>(isBoundary.pos(i,j));
//
//			if(phi(pos) >= zero || cfd.isKinematicBoundary.get_nearest(pos))	isBoundary(i,j) = true;
//			else																isBoundary(i,j) = false;
//		}		
//	}
//	
//	const real h = isBoundary.dx;
//	const double hInv_sq = SQ(1.0/isBoundary.dx);	
//	const double one = 1.0;
//	
//	// build Poisson matrix
//	poisson.zero();
//
//	for ( index_type i = 1; i < poisson.nX-1; ++i ) {
//		Vector2<real> pos, _pos;
//		for ( index_type j = 1; j < poisson.nY-1; ++j ) {
//
//			pos = Vector2<real>(isBoundary.pos(i,j));
//			if( isBoundary(i,j) == false )	{
//
//				_pos = pos - Vector2<real>(h,0); aabb.constrainedBdry(_pos);
//				if( !cfd.isKinematicBoundary.get_nearest(_pos) )	
//					poisson(i,j,0) += one;
//
//				_pos = pos + Vector2<real>(h,0); aabb.constrainedBdry(_pos);
//				if( !cfd.isKinematicBoundary.get_nearest(_pos) )	{
//					poisson(i,j,0) += one;
//					if(phi(_pos) < zero)
//						poisson(i,j,1) -= one;
//				}
//
//				_pos = pos - Vector2<real>(0,h); aabb.constrainedBdry(_pos);
//				if( !cfd.isKinematicBoundary.get_nearest(_pos) )	
//					poisson(i,j,0) += one;
//
//				_pos = pos + Vector2<real>(0,h); aabb.constrainedBdry(_pos);
//				if( !cfd.isKinematicBoundary.get_nearest(_pos) )	{
//					poisson(i,j,0) += one;
//					if(phi(_pos) < zero)
//						poisson(i,j,2) -= one;
//				}
//			}
//		}
//	}	
//}

//void CFDProjection_FreeSurface::buildPoissonMatrix(Grid3_PoissonMatrix<double>& matrix, const Grid3r& phi, const CFDsolver3D& cfd) {
//
//	const real				zero = real(0.0);
//
//	Grid3x4d&				poisson = matrix.poisson;
//	Grid3<bool,double>&		isBoundary = matrix.isBoundary;
//
//	const AABB<real>		aabb = phi.getAABB();
//
//	isBoundary.set_boundary(true);
//
//#ifdef USE_OPENMP
//	int my_end = int(isBoundary.nX-1);
//	int i;
//
//#pragma omp parallel for
//#else
//	index_type my_end = isBoundary.nX-1;
//	index_type i;
//#endif
//	for( i = 1; i < my_end; ++i ) {
//		Vector3<real> pos;
//		for( index_type j = 1; j < isBoundary.nY-1; ++j ) {
//			for( index_type k = 1; k < isBoundary.nZ-1 ; ++k ) {
//				pos = Vector3<real>(isBoundary.pos(i,j,k));
//
//				if(phi(pos) >= zero || cfd.isKinematicBoundary.get_nearest(pos))	isBoundary(i,j,k) = true;
//				else																isBoundary(i,j,k) = false;
//			}
//		}		
//	}
//
//	const real h = isBoundary.dx;
//	const double hInv_sq = SQ(1.0/isBoundary.dx);	
//	const double one = 1.0;
//
//	// build Poisson matrix
//	poisson.zero();
//
//#ifdef USE_OPENMP
//	my_end = int(poisson.nX-1);	
//
//#pragma omp parallel for
//#else
//	my_end = poisson.nX-1;	
//#endif
//	for( i = 1; i < my_end; ++i ) {
//
//		Vector3<real> pos, _pos;
//
//		for ( index_type j = 1; j < poisson.nY-1; ++j ) {
//			for( index_type k = 1; k < poisson.nZ-1; ++k )	{
//
//				pos = Vector3<real>(isBoundary.pos(i,j,k));
//
//				if( isBoundary(i,j,k) == false )	{
//
//					_pos = pos - Vector3<real>(h,0,0); aabb.constrainedBdry(_pos);
//					if( !cfd.isKinematicBoundary.get_nearest(_pos) )							
//						poisson(i,j,k,0) += one;
//
//					_pos = pos + Vector3<real>(h,0,0); aabb.constrainedBdry(_pos);
//					if( !cfd.isKinematicBoundary.get_nearest(_pos) )	{						
//						poisson(i,j,k,0) += one;
//						if(phi(_pos) < zero)
//							poisson(i,j,k,1) -= one;
//					}
//
//					_pos = pos - Vector3<real>(0,h,0); aabb.constrainedBdry(_pos);
//					if( !cfd.isKinematicBoundary.get_nearest(_pos) )							
//						poisson(i,j,k,0) += one;
//
//					_pos = pos + Vector3<real>(0,h,0); aabb.constrainedBdry(_pos);
//					if( !cfd.isKinematicBoundary.get_nearest(_pos) )	{						
//						poisson(i,j,k,0) += one;
//						if(phi(_pos) < zero)
//							poisson(i,j,k,2) -= one;
//					}
//
//					_pos = pos - Vector3<real>(0,0,h); aabb.constrainedBdry(_pos);
//					if( !cfd.isKinematicBoundary.get_nearest(_pos) )							
//						poisson(i,j,k,0) += one;
//
//					_pos = pos + Vector3<real>(0,0,h); aabb.constrainedBdry(_pos);
//					if( !cfd.isKinematicBoundary.get_nearest(_pos) )	{						
//						poisson(i,j,k,0) += one;
//						if(phi(_pos) < zero)
//							poisson(i,j,k,3) -= one;
//					}
//				}
//			}
//		}
//	}
//}
