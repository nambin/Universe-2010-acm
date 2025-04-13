#include <stdafx.h>
#include <Nature/Fluid/CFDsolver/CFDProjection.h>

using namespace std;
using namespace geometryspace;
using namespace fluidspace;

void CFDProjection::init_PoissonSolver_CG(const CFDsolver2D& cfd) { 
	if(sol2D)
		delete sol2D;
	
	sol2D = new Grid2_PoissonSolver_CG<double>(cfd.nX+2, cfd.nY+2, cfd.h, Vector2<double>(-0.5,-0.5));			
}
void CFDProjection::init_PoissonSolver_CG(const CFDsolver3D& cfd) { 
	if(sol3D)
		delete sol3D;

	sol3D = new Grid3_PoissonSolver_CG<double>(cfd.nX+2, cfd.nY+2, cfd.nZ+2, cfd.h, Vector3<double>(-0.5,-0.5,-0.5));			
}

void CFDProjection::setPreconditioner_ICCG() { 
	if(sol2D != NULL && basis::TypeChecker<Grid2_PoissonSolver_CG<double>>()(sol2D)) static_cast<Grid2_PoissonSolver_CG<double>*>(sol2D)->setPreconditioner_ICCG();
	if(sol3D != NULL && basis::TypeChecker<Grid3_PoissonSolver_CG<double>>()(sol3D)) static_cast<Grid3_PoissonSolver_CG<double>*>(sol3D)->setPreconditioner_ICCG();

	log_info << "CG PoissonSolver _ ICCG" << endl;
}

void CFDProjection::setPreconditioner_BlockRedBlackICCG(index_type nbX, index_type nbY)	{ 
	if(sol2D != NULL && basis::TypeChecker<Grid2_PoissonSolver_CG<double>>()(sol2D)) 
		static_cast<Grid2_PoissonSolver_CG<double>*>(sol2D)->setPreconditioner_BlockRedBlackICCG(nbX, nbY);

	log_info << "CG PoissonSolver _ BlockRedBlackICCG " << "(" << nbX << ", " << nbY << ")" << endl;
}
void CFDProjection::setPreconditioner_BlockRedBlackICCG(index_type nbX, index_type nbY, index_type nbZ)	{ 
	if(sol3D != NULL && basis::TypeChecker<Grid3_PoissonSolver_CG<double>>()(sol3D)) 
		static_cast<Grid3_PoissonSolver_CG<double>*>(sol3D)->setPreconditioner_BlockRedBlackICCG(nbX, nbY, nbZ);	

	log_info << "CG PoissonSolver _ BlockRedBlackICCG " << "(" << nbX << ", " << nbY << ", " << nbZ << ")" << endl;
}

void CFDProjection::init_PoissonSolver_GMG(const CFDsolver2D& cfd, unsigned int num_level) { 
	if(sol2D)
		delete sol2D;

	sol2D = new Grid2_PoissonSolver_GMG<double>(cfd.nX+2, cfd.nY+2, cfd.h, Vector2<double>(-0.5,-0.5), num_level);			
}
void CFDProjection::init_PoissonSolver_GMG(const CFDsolver3D& cfd, unsigned int num_level) { 
	if(sol3D)
		delete sol3D;

	sol3D = new Grid3_PoissonSolver_GMG<double>(cfd.nX+2, cfd.nY+2, cfd.nZ+2, cfd.h, Vector3<double>(-0.5,-0.5,-0.5), num_level);			
}

void CFDProjection::set_iteration_fine_grid_GMG(size_t iter) {
	if(sol2D != NULL && basis::TypeChecker<Grid2_PoissonSolver_GMG<double>>()(sol2D)) static_cast<Grid2_PoissonSolver_GMG<double>*>(sol2D)->set_iteration_number_fine_grid(iter);
	if(sol3D != NULL && basis::TypeChecker<Grid3_PoissonSolver_GMG<double>>()(sol3D)) static_cast<Grid3_PoissonSolver_GMG<double>*>(sol3D)->set_iteration_number_fine_grid(iter);
}

void CFDProjection::set_iteration_number_GMG(size_t pre, size_t post) {
	if(sol2D != NULL && basis::TypeChecker<Grid2_PoissonSolver_GMG<double>>()(sol2D)) static_cast<Grid2_PoissonSolver_GMG<double>*>(sol2D)->set_iteration_number(pre, post);
	if(sol3D != NULL && basis::TypeChecker<Grid3_PoissonSolver_GMG<double>>()(sol3D)) static_cast<Grid3_PoissonSolver_GMG<double>*>(sol3D)->set_iteration_number(pre, post);
}

void CFDProjection::set_parallel_threshold_GMG(index_type thres) {
	if(sol2D != NULL && basis::TypeChecker<Grid2_PoissonSolver_GMG<double>>()(sol2D)) static_cast<Grid2_PoissonSolver_GMG<double>*>(sol2D)->set_parallel_threshold(thres);
	if(sol3D != NULL && basis::TypeChecker<Grid3_PoissonSolver_GMG<double>>()(sol3D)) static_cast<Grid3_PoissonSolver_GMG<double>*>(sol3D)->set_parallel_threshold(thres);
}

void CFDProjection::setComputingOption(const ParallelComputable& rhs) {
	ParallelComputable::setComputingOption(rhs);
	if(sol2D) sol2D->setComputingOption(rhs);
	if(sol3D) sol3D->setComputingOption(rhs);
}

void CFDProjection::buildDivergence(Grid2d &divergence, const CFDsolver2D& cfd) {
	
	const Grid2b& isKinematicBoundary = cfd.isKinematicBoundary;		
	const Grid2r& u = cfd.u;
	const Grid2r& v = cfd.v;	
	
	divergence.zero();
	for( index_type i=1; i != divergence.nX-1; ++i ) {
		for( index_type j=1; j != divergence.nY-1; ++j ) {
			if(!isKinematicBoundary(i,j))	{
//				divergence(i,j) = ( ( u(cfd.iS2U(i,j,1)) - u(cfd.iS2U(i,j,0)) ) 
//								  + ( v(cfd.iS2V(i,j,1)) - v(cfd.iS2V(i,j,0)) ) ) * cfd.hInv;

				divergence(i,j) = ( ( u(cfd.iS2U(i,j,1)) - u(cfd.iS2U(i,j,0)) ) 
								  + ( v(cfd.iS2V(i,j,1)) - v(cfd.iS2V(i,j,0)) ) ) ;	
			}					
		}
	}
}

void CFDProjection::buildDivergence(Grid3d &divergence, const CFDsolver3D& cfd) {

	const Grid3b& isKinematicBoundary = cfd.isKinematicBoundary;		
	const Grid3r& u = cfd.u, & v = cfd.v, & w = cfd.w;	

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
		for( index_type j = 1; j < divergence.nY-1; ++j ) {
			for( index_type k = 1; k < divergence.nZ-1; ++k ) {
				if(!isKinematicBoundary(i,j,k)) {
					//divergence(i,j,k) = ( ( u(cfd.iS2U(i,j,k,1)) - u(cfd.iS2U(i,j,k,0)) )
					//					+ ( v(cfd.iS2V(i,j,k,1)) - v(cfd.iS2V(i,j,k,0)) )
					//					+ ( w(cfd.iS2W(i,j,k,1)) - w(cfd.iS2W(i,j,k,0)) ) )  * cfd.hInv;

					divergence(i,j,k) = ( ( u(cfd.iS2U(i,j,k,1)) - u(cfd.iS2U(i,j,k,0)) )
										+ ( v(cfd.iS2V(i,j,k,1)) - v(cfd.iS2V(i,j,k,0)) )
										+ ( w(cfd.iS2W(i,j,k,1)) - w(cfd.iS2W(i,j,k,0)) ) )  ;
				}					
			}
		}
	}	
}

void CFDProjection::addPressureGradient(CFDsolver2D& cfd, const Grid2d& pressure) {

	const Grid2b&		isKinematicBoundary = cfd.isKinematicBoundary;
	Grid2r&				u = cfd.u;
	Grid2r&				v = cfd.v;	

	for( index_type i=1; i < u.nX-1; ++i ) {
		for( index_type j=1; j < u.nY-1; ++j ) {					
			if( !isKinematicBoundary(cfd.iU2S(i,j,1)) && !isKinematicBoundary(cfd.iU2S(i,j,0)) ) {
				//u(i,j) += ( pressure(cfd.iU2S(i,j,1)) - pressure(cfd.iU2S(i,j,0)) ) * cfd.hInv;
				u(i,j) += ( pressure(cfd.iU2S(i,j,1)) - pressure(cfd.iU2S(i,j,0)) ) ;
			}
		}
	}

	for( index_type i=1; i < v.nX-1; ++i ) {
		for( index_type j=1; j < v.nY-1; ++j ) {					
			if( !isKinematicBoundary(cfd.iV2S(i,j,1)) && !isKinematicBoundary(cfd.iV2S(i,j,0)) ) {
				//v(i,j) += ( pressure(cfd.iV2S(i,j,1)) - pressure(cfd.iV2S(i,j,0)) ) * cfd.hInv;
				v(i,j) += ( pressure(cfd.iV2S(i,j,1)) - pressure(cfd.iV2S(i,j,0)) ) ;
			}
		}
	}
}

void CFDProjection::addPressureGradient(CFDsolver3D& cfd, const Grid3d& pressure) {

#ifdef USE_OPENMP	
	int my_end = int(cfd.u.nX-1);
	int i;
#pragma omp parallel for
#else	
	index_type my_end = cfd.u.nX-1;
	index_type i;
#endif
	for( i = 1; i < my_end; ++i ) {	
		for ( index_type j = 1; j < cfd.u.nY-1; ++j ) {
			for( index_type k = 1; k < cfd.u.nZ-1; ++k )	{
				if( !cfd.isKinematicBoundary[cfd.iU2S(i,j,k,1)] && !cfd.isKinematicBoundary[cfd.iU2S(i,j,k,0)] )
				{				
					//cfd.u(i,j,k) += ( pressure[cfd.iU2S(i,j,k,1)] - pressure[cfd.iU2S(i,j,k,0)] ) * cfd.hInv;
					cfd.u(i,j,k) += ( pressure[cfd.iU2S(i,j,k,1)] - pressure[cfd.iU2S(i,j,k,0)] ) ;
				}
			}
		}
	}

#ifdef USE_OPENMP	
	my_end = int(cfd.v.nX-1);
#pragma omp parallel for
#else	
	my_end = cfd.v.nX-1;
#endif
	for( i = 1; i < my_end; ++i ) {
		for ( index_type j = 1; j < cfd.v.nY-1; ++j ) {
			for( index_type k = 1; k < cfd.v.nZ-1; ++k ) {
				if( !cfd.isKinematicBoundary[cfd.iV2S(i,j,k,1)] && !cfd.isKinematicBoundary[cfd.iV2S(i,j,k,0)] )	
				{
					//cfd.v(i,j,k) += ( pressure[cfd.iV2S(i,j,k,1)] - pressure[cfd.iV2S(i,j,k,0)] ) * cfd.hInv;
					cfd.v(i,j,k) += ( pressure[cfd.iV2S(i,j,k,1)] - pressure[cfd.iV2S(i,j,k,0)] ) ;
				}
			}
		}
	}

#ifdef USE_OPENMP	
	my_end = int(cfd.w.nX-1);
#pragma omp parallel for
#else	
	my_end = cfd.w.nX-1;
#endif
	for( i = 1; i < my_end; ++i ) {	
		for ( index_type j = 1; j < cfd.w.nY-1; ++j ) {
			for( index_type k = 1; k < cfd.w.nZ-1; ++k ) {
				if( !cfd.isKinematicBoundary[cfd.iW2S(i,j,k,1)] && !cfd.isKinematicBoundary[cfd.iW2S(i,j,k,0)] )	
				{
					//cfd.w(i,j,k) += ( pressure[cfd.iW2S(i,j,k,1)] - pressure[cfd.iW2S(i,j,k,0)] ) * cfd.hInv;
					cfd.w(i,j,k) += ( pressure[cfd.iW2S(i,j,k,1)] - pressure[cfd.iW2S(i,j,k,0)] ) ;
				}
			}
		}
	}
}

void CFDProjection::buildPoissonMatrix(Grid2x3d& poisson, Grid2<bool,double>& isbdry, const Grid2d& u_part, const Grid2d& v_part) {

	assert(poisson.nX == u_part.nX + 1);	assert(poisson.nY == u_part.nY);
	assert(poisson.nX == v_part.nX);		assert(poisson.nY == v_part.nY + 1);
		
	poisson.zero();	
	for( index_type i=1; i < poisson.nX-1 ; ++i ) {
		for( index_type j=1; j < poisson.nY-1 ; ++j ) {			
			poisson(i,j,0) = u_part(i-1,j) + u_part(i,j) + v_part(i,j-1) + v_part(i,j);
			poisson(i,j,1) = - u_part(i,j);
			poisson(i,j,2) = - v_part(i,j);			
		}
	}	

	for( index_type i=0; i < poisson.nX ; ++i ) {
		for( index_type j=0; j < poisson.nY ; ++j ) {
			if(poisson(i,j,0) == 0) isbdry(i,j) = true;
			else					isbdry(i,j) = false;
		}
	}
}

void CFDProjection::buildPoissonMatrix(Grid3x4d& poisson, Grid3<bool,double>& isbdry, const Grid3d& u_part, const Grid3d& v_part, const Grid3d& w_part) {
	
	assert(poisson.nX == u_part.nX + 1);	assert(poisson.nY == u_part.nY);		assert(poisson.nZ == u_part.nZ);
	assert(poisson.nX == v_part.nX);		assert(poisson.nY == v_part.nY + 1);	assert(poisson.nZ == v_part.nZ);
	assert(poisson.nX == w_part.nX);		assert(poisson.nY == w_part.nY);		assert(poisson.nZ == w_part.nZ + 1);

	poisson.zero();	

#ifdef USE_OPENMP
	int my_end = int(poisson.nX - 1);
	int i;

#pragma omp parallel for
#else
	index_type my_end = poisson.nX - 1;
	index_type i;
#endif
	for( i = 1; i < my_end; ++i ) {	
		for( index_type j = 1; j < poisson.nY - 1; ++j ) {
			for( index_type k = 1; k < poisson.nZ - 1; ++k ) {
				poisson(i,j,k,0) = u_part(i-1,j,k) + u_part(i,j,k) + v_part(i,j-1,k) + v_part(i,j,k) + w_part(i,j,k-1) + w_part(i,j,k);
				poisson(i,j,k,1) = - u_part(i,j,k);
				poisson(i,j,k,2) = - v_part(i,j,k);
				poisson(i,j,k,3) = - w_part(i,j,k);								
			}
		}
	}

#ifdef USE_OPENMP
	my_end = int(poisson.nX);

#pragma omp parallel for
#else
	my_end = poisson.nX;
#endif
	for( i = 0; i < my_end; ++i ) {	
		for( index_type j = 0; j < poisson.nY ; ++j ) {
			for( index_type k = 0; k < poisson.nZ ; ++k ) {
				if(poisson(i,j,k,0) == 0)	isbdry(i,j,k) = true;
				else						isbdry(i,j,k) = false;
			}
		}
	}	
}

void CFDProjection::buildPoissonMatrix(Grid2x3d& poisson, Grid2<bool,double>& isbdry, const Grid2d& u_part, const Grid2d& v_part, const Grid2d& s_part) {
	assert(poisson.nX == u_part.nX + 1);	assert(poisson.nY == u_part.nY);
	assert(poisson.nX == v_part.nX);		assert(poisson.nY == v_part.nY + 1);
	assert(poisson.nX == s_part.nX);		assert(poisson.nY == s_part.nY);

	poisson.zero();	
	for( index_type i=1; i < poisson.nX-1 ; ++i ) {
		for( index_type j=1; j < poisson.nY-1 ; ++j ) {			
			poisson(i,j,0) = u_part(i-1,j) + u_part(i,j) + v_part(i,j-1) + v_part(i,j) + s_part(i,j);
			poisson(i,j,1) = - u_part(i,j);
			poisson(i,j,2) = - v_part(i,j);			
		}
	}	

	for( index_type i=0; i < poisson.nX ; ++i ) {
		for( index_type j=0; j < poisson.nY ; ++j ) {
			if(poisson(i,j,0) == 0) isbdry(i,j) = true;
			else					isbdry(i,j) = false;
		}
	}
}

void CFDProjection::buildPoissonMatrix(Grid3x4d& poisson, Grid3<bool,double>& isbdry, const Grid3d& u_part, const Grid3d& v_part, const Grid3d& w_part, const Grid3d& s_part) {
	assert(poisson.nX == u_part.nX + 1);	assert(poisson.nY == u_part.nY);		assert(poisson.nZ == u_part.nZ);
	assert(poisson.nX == v_part.nX);		assert(poisson.nY == v_part.nY + 1);	assert(poisson.nZ == v_part.nZ);
	assert(poisson.nX == w_part.nX);		assert(poisson.nY == w_part.nY);		assert(poisson.nZ == w_part.nZ + 1);
	assert(poisson.nX == s_part.nX);		assert(poisson.nY == s_part.nY);		assert(poisson.nZ == s_part.nZ);

	poisson.zero();	

#ifdef USE_OPENMP
	int my_end = int(poisson.nX - 1);
	int i;

#pragma omp parallel for
#else
	index_type my_end = poisson.nX - 1;
	index_type i;
#endif
	for( i = 1; i < my_end; ++i ) {		
		for( index_type j=1; j < poisson.nY-1 ; ++j ) {			
			for( index_type k=1; k < poisson.nZ-1 ; ++k ) {		
				poisson(i,j,k,0) = u_part(i-1,j,k) + u_part(i,j,k) + v_part(i,j-1,k) + v_part(i,j,k) + w_part(i,j,k-1) + w_part(i,j,k) + s_part(i,j,k);
				poisson(i,j,k,1) = - u_part(i,j,k);
				poisson(i,j,k,2) = - v_part(i,j,k);			
				poisson(i,j,k,3) = - w_part(i,j,k);
			}			
		}
	}	

#ifdef USE_OPENMP
	my_end = int(poisson.nX);

#pragma omp parallel for
#else
	my_end = poisson.nX;
#endif
	for( i = 0; i < my_end; ++i ) {		
		for( index_type j=0; j < poisson.nY ; ++j ) {
			for( index_type k=0; k < poisson.nZ ; ++k ) {
				if(poisson(i,j,k,0) == 0)	isbdry(i,j,k) = true;
				else						isbdry(i,j,k) = false;
			}			
		}
	}
}

void CFDProjection::coarsening_PoissonMatrix(Grid2d& u_coarse, Grid2d& v_coarse, const Grid2d& u_fine, const Grid2d& v_fine) {
	const index_type nX = u_fine.nX - 1;
	const index_type nY = u_fine.nY - 2;

	const double dx = u_fine.dx;
	const double dy = u_fine.dy;

	if(nX % 2 != 0 || nY % 2 != 0)
		throwError("Error in CFDProjection::coarsening_PoissonMatrix()");

	u_coarse.init( nX / 2 + 1, nY / 2 + 2, dx*2, dy*2, u_fine.offset );	u_coarse.zero();
	v_coarse.init( nX / 2 + 2, nY / 2 + 1, dx*2, dy*2, v_fine.offset );	v_coarse.zero();

	const double eighth = double(1.0 / 8.0);

	for(index_type i=1;i<u_coarse.nX-1;++i)
		for(index_type j=1;j<u_coarse.nY-1;++j) 
			u_coarse(i,j) = (u_fine(i*2, j*2-1) + u_fine(i*2, j*2)) * eighth;

	for(index_type i=1;i<v_coarse.nX-1;++i)
		for(index_type j=1;j<v_coarse.nY-1;++j) 
			v_coarse(i,j) = (v_fine(i*2-1, j*2) + v_fine(i*2, j*2)) * eighth;
}

void CFDProjection::coarsening_PoissonMatrix(Grid2d& u_coarse, Grid2d& v_coarse, Grid2d& s_coarse, const Grid2d& u_fine, const Grid2d& v_fine, const Grid2d& s_fine) {
	coarsening_PoissonMatrix(u_coarse, v_coarse, u_fine, v_fine);

	const index_type nX = s_fine.nX - 2;
	const index_type nY = s_fine.nY - 2;

	const double dx = s_fine.dx;
	const double dy = s_fine.dy;

	if(nX % 2 != 0 || nY % 2 != 0)
		throwError("Error in CFDProjection::coarsening_PoissonMatrix()");

	s_coarse.init( nX / 2 + 2, nY / 2 + 2, dx*2, dy*2, s_fine.offset );	s_coarse.zero();

	const double divisor = double(1.0 / 4.0);

	for(index_type i=1;i<s_coarse.nX-1;++i)
		for(index_type j=1;j<s_coarse.nY-1;++j) 
			s_coarse(i,j) = ( s_fine(i*2-1, j*2-1) + s_fine(i*2, j*2-1) 
							+ s_fine(i*2-1, j*2  ) + s_fine(i*2, j*2  )) * divisor;
}

void CFDProjection::coarsening_PoissonMatrix(Grid3d& u_coarse, Grid3d& v_coarse, Grid3d& w_coarse, const Grid3d& u_fine, const Grid3d& v_fine, const Grid3d& w_fine) {
	const index_type nX = u_fine.nX - 1;
	const index_type nY = u_fine.nY - 2;
	const index_type nZ = u_fine.nZ - 2;

	const double dx = u_fine.dx;
	const double dy = u_fine.dy;
	const double dz = u_fine.dz;

	if(nX % 2 != 0 || nY % 2 != 0 || nZ % 2 != 0)
		throwError("Error in CFDProjection::coarsening_PoissonMatrix()");

	u_coarse.init( nX / 2 + 1, nY / 2 + 2, nZ / 2 + 2, dx*2, dy*2, dz*2, u_fine.offset );	u_coarse.zero();
	v_coarse.init( nX / 2 + 2, nY / 2 + 1, nZ / 2 + 2, dx*2, dy*2, dz*2, v_fine.offset );	v_coarse.zero();
	w_coarse.init( nX / 2 + 2, nY / 2 + 2, nZ / 2 + 1, dx*2, dy*2, dz*2, w_fine.offset );	w_coarse.zero();

	const double sixteenth = double(1.0 / 16.0);

	for(index_type i=1;i<u_coarse.nX-1;++i)
		for(index_type j=1;j<u_coarse.nY-1;++j) 
			for(index_type k=1;k<u_coarse.nZ-1;++k) 
				u_coarse(i,j,k) = (	u_fine(i*2, j*2-1, k*2-1) + u_fine(i*2, j*2, k*2-1)
								  + u_fine(i*2, j*2-1, k*2  ) + u_fine(i*2, j*2, k*2  ) ) * sixteenth;

	for(index_type i=1;i<v_coarse.nX-1;++i)
		for(index_type j=1;j<v_coarse.nY-1;++j) 
			for(index_type k=1;k<v_coarse.nZ-1;++k) 
				v_coarse(i,j,k) = (	v_fine(i*2-1, j*2, k*2-1) + v_fine(i*2, j*2, k*2-1)
								  + v_fine(i*2-1, j*2, k*2  ) + v_fine(i*2, j*2, k*2  )	) * sixteenth;

	for(index_type i=1;i<w_coarse.nX-1;++i)
		for(index_type j=1;j<w_coarse.nY-1;++j) 
			for(index_type k=1;k<w_coarse.nZ-1;++k) 
				w_coarse(i,j,k) = (	w_fine(i*2-1, j*2-1, k*2  ) + w_fine(i*2, j*2-1, k*2  )
								  + w_fine(i*2-1, j*2  , k*2  ) + w_fine(i*2, j*2  , k*2  )	) * sixteenth;
}

void CFDProjection::coarsening_PoissonMatrix(Grid3d& u_coarse, Grid3d& v_coarse, Grid3d& w_coarse, Grid3d& s_coarse, const Grid3d& u_fine, const Grid3d& v_fine, const Grid3d& w_fine, const Grid3d& s_fine) {
	coarsening_PoissonMatrix(u_coarse, v_coarse, w_coarse, u_fine, v_fine, w_fine);

	const index_type nX = s_fine.nX - 2;
	const index_type nY = s_fine.nY - 2;
	const index_type nZ = s_fine.nZ - 2;

	const double dx = s_fine.dx;
	const double dy = s_fine.dy;
	const double dz = s_fine.dz;

	if(nX % 2 != 0 || nY % 2 != 0)
		throwError("Error in CFDProjection::coarsening_PoissonMatrix()");

	s_coarse.init( nX / 2 + 2, nY / 2 + 2, nZ / 2 + 2, dx*2, dy*2, dz*2, s_fine.offset );	s_coarse.zero();

	const double divisor = double(1.0 / 8.0);

	for(index_type i=1;i<s_coarse.nX-1;++i)
		for(index_type j=1;j<s_coarse.nY-1;++j) 
			for(index_type k=1;k<s_coarse.nZ-1;++k) 
				s_coarse(i,j,k) = ( s_fine(i*2-1, j*2-1, k*2-1) + s_fine(i*2, j*2-1, k*2-1) + s_fine(i*2-1, j*2, k*2-1) + s_fine(i*2, j*2, k*2-1) 
								  + s_fine(i*2-1, j*2-1, k*2  ) + s_fine(i*2, j*2-1, k*2  ) + s_fine(i*2-1, j*2, k*2  ) + s_fine(i*2, j*2, k*2  ) ) * divisor;
}

void CFDProjection::buildPoissonMatrix(Grid2_PoissonSolver_GMG<double>& gmg, const Grid2d& u_part, const Grid2d& v_part) {
	buildPoissonMatrix(gmg.poissonMatrix.poisson, gmg.poissonMatrix.isBoundary, u_part, v_part);

	Grid2d* u_fine = new Grid2d(), *v_fine = new Grid2d();
	Grid2d* u_coarse = new Grid2d(), *v_coarse = new Grid2d();

	for(size_t lv=1;lv<gmg.solvers.size();++lv) {
		if(lv == 1) coarsening_PoissonMatrix(*u_coarse, *v_coarse, u_part, v_part);
		else		coarsening_PoissonMatrix(*u_coarse, *v_coarse, *u_fine, *v_fine);

		buildPoissonMatrix(gmg.solvers[lv]->poissonMatrix.poisson, gmg.solvers[lv]->poissonMatrix.isBoundary, *u_coarse, *v_coarse);

		swap(u_fine, u_coarse);
		swap(v_fine, v_coarse);
	}

	delete u_fine;	delete u_coarse;	
	delete v_fine;  delete v_coarse;
}

void CFDProjection::buildPoissonMatrix(Grid3_PoissonSolver_GMG<double>& gmg, const Grid3d& u_part, const Grid3d& v_part, const Grid3d& w_part) {
	buildPoissonMatrix(gmg.poissonMatrix.poisson, gmg.poissonMatrix.isBoundary, u_part, v_part, w_part);

	Grid3d* u_fine = new Grid3d(), *v_fine = new Grid3d(), *w_fine = new Grid3d();
	Grid3d* u_coarse = new Grid3d(), *v_coarse = new Grid3d(), *w_coarse = new Grid3d();

	for(size_t lv=1;lv<gmg.solvers.size();++lv) {
		if(lv == 1) coarsening_PoissonMatrix(*u_coarse, *v_coarse, *w_coarse, u_part, v_part, w_part);
		else		coarsening_PoissonMatrix(*u_coarse, *v_coarse, *w_coarse, *u_fine, *v_fine, *w_fine);

		buildPoissonMatrix(gmg.solvers[lv]->poissonMatrix.poisson, gmg.solvers[lv]->poissonMatrix.isBoundary, *u_coarse, *v_coarse, *w_coarse);

		swap(u_fine, u_coarse);
		swap(v_fine, v_coarse);
		swap(w_fine, w_coarse);
	}

	delete u_fine;	delete u_coarse;	
	delete v_fine;  delete v_coarse;
	delete w_fine;  delete w_coarse;
}

void CFDProjection::buildPoissonMatrix(Grid2_PoissonSolver_GMG<double>& gmg, const Grid2d& u_part, const Grid2d& v_part, const Grid2d& s_part) {
	buildPoissonMatrix(gmg.poissonMatrix.poisson, gmg.poissonMatrix.isBoundary, u_part, v_part, s_part);

	Grid2d* u_fine = new Grid2d(), *v_fine = new Grid2d(), *s_fine = new Grid2d();
	Grid2d* u_coarse = new Grid2d(), *v_coarse = new Grid2d(), *s_coarse = new Grid2d();	

	for(size_t lv=1;lv<gmg.solvers.size();++lv) {
		if(lv == 1) coarsening_PoissonMatrix(*u_coarse, *v_coarse, *s_coarse, u_part, v_part, s_part);
		else		coarsening_PoissonMatrix(*u_coarse, *v_coarse, *s_coarse, *u_fine, *v_fine, *s_fine);

		buildPoissonMatrix(gmg.solvers[lv]->poissonMatrix.poisson, gmg.solvers[lv]->poissonMatrix.isBoundary, *u_coarse, *v_coarse, *s_coarse);

		swap(u_fine, u_coarse);
		swap(v_fine, v_coarse);
		swap(s_fine, s_coarse);
	}

	delete u_fine;	delete u_coarse;	
	delete v_fine;  delete v_coarse;
	delete s_fine;  delete s_coarse;
}

void CFDProjection::buildPoissonMatrix(Grid3_PoissonSolver_GMG<double>& gmg, const Grid3d& u_part, const Grid3d& v_part, const Grid3d& w_part, const Grid3d& s_part) {
	buildPoissonMatrix(gmg.poissonMatrix.poisson, gmg.poissonMatrix.isBoundary, u_part, v_part, w_part, s_part);

	Grid3d* u_fine = new Grid3d(), *v_fine = new Grid3d(), *w_fine = new Grid3d(), *s_fine = new Grid3d();
	Grid3d* u_coarse = new Grid3d(), *v_coarse = new Grid3d(), *w_coarse = new Grid3d(), *s_coarse = new Grid3d();	

	for(size_t lv=1;lv<gmg.solvers.size();++lv) {
		if(lv == 1) coarsening_PoissonMatrix(*u_coarse, *v_coarse, *w_coarse, *s_coarse, u_part, v_part, w_part, s_part);
		else		coarsening_PoissonMatrix(*u_coarse, *v_coarse, *w_coarse, *s_coarse, *u_fine, *v_fine, *w_fine, *s_fine);

		buildPoissonMatrix(gmg.solvers[lv]->poissonMatrix.poisson, gmg.solvers[lv]->poissonMatrix.isBoundary, *u_coarse, *v_coarse, *w_coarse, *s_coarse);

		swap(u_fine, u_coarse);
		swap(v_fine, v_coarse);
		swap(w_fine, w_coarse);
		swap(s_fine, s_coarse);
	}

	delete u_fine;	delete u_coarse;	
	delete v_fine;  delete v_coarse;
	delete w_fine;  delete w_coarse;
	delete s_fine;  delete s_coarse;
}

void CFDProjection::importBinary(std::ifstream& ifs) {
	if(sol2D) sol2D->x.importBinary(ifs);
	if(sol3D) sol3D->x.importBinary(ifs);
}

void CFDProjection::exportBinary(std::ofstream& ofs) const {
	if(sol2D) sol2D->x.exportBinary(ofs);
	if(sol3D) sol3D->x.exportBinary(ofs);
}