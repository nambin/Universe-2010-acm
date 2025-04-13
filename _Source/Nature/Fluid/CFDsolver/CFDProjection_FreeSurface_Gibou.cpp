#include <stdafx.h>
#include <Nature/Fluid/CFDsolver/CFDProjection.h>

using namespace std;
using namespace fluidspace;

void CFDProjection_FreeSurface_Gibou::buildPoissonMatrix(Grid2_PoissonSolver<double>& solver, const Grid2r& phi, const CFDsolver2D& cfd) {
	if(applySurfaceTension)
		div2D = &solver.b;

	Grid2d u_part, v_part, s_part;
	buildPoissonMatrix(u_part, v_part, s_part, phi, cfd);

	if(TypeChecker<Grid2_PoissonSolver_CG<double>>()(&solver)) 
		CFDProjection::buildPoissonMatrix(solver.poissonMatrix.poisson, solver.poissonMatrix.isBoundary, u_part, v_part, s_part);

	else if(TypeChecker<Grid2_PoissonSolver_GMG<double>>()(&solver)) 
		CFDProjection::buildPoissonMatrix(static_cast<Grid2_PoissonSolver_GMG<double>&>(solver), u_part, v_part, s_part);

	div2D = NULL;
}

void CFDProjection_FreeSurface_Gibou::buildPoissonMatrix(Grid3_PoissonSolver<double>& solver, const Grid3r& phi, const CFDsolver3D& cfd) {
	if(applySurfaceTension)
		div3D = &solver.b;

	Grid3d u_part, v_part, w_part, s_part;
	buildPoissonMatrix(u_part, v_part, w_part, s_part, phi, cfd);

	if(TypeChecker<Grid3_PoissonSolver_CG<double>>()(&solver)) 
		CFDProjection::buildPoissonMatrix(solver.poissonMatrix.poisson, solver.poissonMatrix.isBoundary, u_part, v_part, w_part, s_part);

	else if(TypeChecker<Grid3_PoissonSolver_GMG<double>>()(&solver)) 
		CFDProjection::buildPoissonMatrix(static_cast<Grid3_PoissonSolver_GMG<double>&>(solver), u_part, v_part, w_part, s_part);

	div3D = NULL;
}

void CFDProjection_FreeSurface_Gibou::buildDivergence(Grid2d& divergence, const Grid2r& phi, const CFDsolver2D& cfd) {

	const real zero = 0.0;

	if(applySurfaceTension) {

		Grid2r& u = const_cast<Grid2r&>(cfd.u);
		Grid2r& v = const_cast<Grid2r&>(cfd.v);

		curvature2D.resize(cfd.isKinematicBoundary);
		tracking->getCurvature(curvature2D);
	
		const real tension_factor = sigma * dt * cfd.hInv;

		for( index_type i=1; i < cfd.u.nX-1 ; ++i ) {
			index_type rt, lt;
			double theta, inv_theta;

			for( index_type j=1; j < cfd.u.nY-1 ; ++j ) {					
				rt = cfd.iU2S(i,j,1);
				lt = cfd.iU2S(i,j,0);

				if( !cfd.isKinematicBoundary(lt) && !cfd.isKinematicBoundary(rt) ) {	

					if( phi(lt) > zero && phi(rt) < zero ) {	
						theta = phi(rt) / (phi(rt) - phi(lt));
						inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));
						
						const double curvature = lerp(curvature2D(lt), curvature2D(rt), theta);
						const double factor = tension_factor * inv_theta * curvature;

						u(i,j) += factor;						
					}
					else if( phi(lt) < zero && phi(rt) > zero ) {
						theta = phi(lt) / (phi(lt) - phi(rt));
						inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

						const double curvature = lerp(curvature2D(rt), curvature2D(lt), theta);
						const double factor = tension_factor * inv_theta * curvature;

						u(i,j) -= factor;						
					}
				}
			}
		}

		for( index_type i=1; i < cfd.v.nX-1 ; ++i ) {
			index_type rt, lt;
			double theta, inv_theta;

			for( index_type j=1; j < cfd.v.nY-1 ; ++j ) {					
				rt = cfd.iV2S(i,j,1);
				lt = cfd.iV2S(i,j,0);

				if( !cfd.isKinematicBoundary(lt) && !cfd.isKinematicBoundary(rt) ) {
					
					if( phi(lt) > zero && phi(rt) < zero ) {
						theta = phi(rt) / (phi(rt) - phi(lt));
						inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

						const double curvature = lerp(curvature2D(lt), curvature2D(rt), theta);
						const double factor = tension_factor * inv_theta * curvature;

						v(i,j) += factor;
					}
					else if( phi(lt) < zero && phi(rt) > zero ) {
						theta = phi(lt) / (phi(lt) - phi(rt));
						inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

						const double curvature = lerp(curvature2D(rt), curvature2D(lt), theta);
						const double factor = tension_factor * inv_theta * curvature;

						v(i,j) -= factor;
					}
				}
			}
		}	
	}	
	
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

void CFDProjection_FreeSurface_Gibou::buildPoissonMatrix(Grid2d& u_part, Grid2d& v_part, Grid2d& s_part, const Grid2r& phi, const CFDsolver2D& cfd) {

	if(u_part.is_empty())	u_part.resize(cfd.u);
	if(v_part.is_empty())	v_part.resize(cfd.v);
	if(s_part.is_empty())	s_part.resize(cfd.isKinematicBoundary);

	const Grid2b&	isKinematicBoundary = cfd.isKinematicBoundary;
	const double	one = 1.0;

	const real zero = real(0.0);

	u_part.zero();
	v_part.zero();
	s_part.zero();

	for ( index_type i = 1; i < isKinematicBoundary.nX-1; ++i ) {
		double theta, inv_theta;

		for ( index_type j = 1; j < isKinematicBoundary.nY-1; ++j ) {
		
			if( phi(i,j) < zero && !isKinematicBoundary(i,j) )	{
				if( !isKinematicBoundary(i-1,j) && phi(i-1,j) >= zero ) {
					theta = phi(i,j) / (phi(i,j) - phi(i-1,j));
					inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

					s_part(i,j) +=  inv_theta;
				}

				if( !isKinematicBoundary(i+1,j) )	{
					if(phi(i+1,j) < zero) {						
						u_part(i,j) += one;
					}
					else {
						theta = phi(i,j) / (phi(i,j) - phi(i+1,j));
						inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

						s_part(i,j) += inv_theta;
					}
				}

				if( !isKinematicBoundary(i,j-1) && phi(i,j-1) >= zero ) {					
					theta = phi(i,j) / (phi(i,j) - phi(i,j-1));
					inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

					s_part(i,j) += inv_theta;
				}

				if( !isKinematicBoundary(i,j+1) )	{					
					if(phi(i,j+1) < zero) {						
						v_part(i,j) += one;
					}
					else {
						theta = phi(i,j) / (phi(i,j) - phi(i,j+1));
						inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

						s_part(i,j) += inv_theta;						
					}
				}
			}
		}
	}
}

void CFDProjection_FreeSurface_Gibou::buildDivergence(Grid3d& divergence, const Grid3r& phi, const CFDsolver3D& cfd) {

	const real zero	= 0.0;	

	if(applySurfaceTension) {

		Grid3r& u = const_cast<Grid3r&>(cfd.u);
		Grid3r& v = const_cast<Grid3r&>(cfd.v);
		Grid3r& w = const_cast<Grid3r&>(cfd.w);

		curvature3D.resize(cfd.isKinematicBoundary);
		tracking->getCurvature(curvature3D);

		const real tension_factor = sigma * dt * cfd.hInv;

#ifdef USE_OPENMP
		int my_end = int(cfd.u.nX-1);
		int i;

#pragma omp parallel for
#else
		index_type my_end = cfd.u.nX-1;
		index_type i;
#endif
		for( i = 1; i < my_end; ++i ) {	
			index_type rt, lt;
			double theta, inv_theta;

			for( index_type j=1; j < cfd.u.nY-1 ; ++j ) {					
				for( index_type k=1; k < cfd.u.nZ-1 ; ++k ) {
					rt = cfd.iU2S(i,j,k,1);
					lt = cfd.iU2S(i,j,k,0);

					if( !cfd.isKinematicBoundary(lt) && !cfd.isKinematicBoundary(rt)) {
						
						if( phi(lt) > zero && phi(rt) < zero ) {
							theta = phi(rt) / (phi(rt) - phi(lt));
							inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

							u(i,j,k) += tension_factor * inv_theta * lerp(curvature3D(lt), curvature3D(rt), theta);
						}
						else if( phi(lt) < zero && phi(rt) > zero ) {
							theta = phi(lt) / (phi(lt) - phi(rt));
							inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

							u(i,j,k) -= tension_factor * inv_theta * lerp(curvature3D(rt), curvature3D(lt), theta);
						}
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
			index_type rt, lt;
			double theta, inv_theta;

			for( index_type j=1; j < cfd.v.nY-1 ; ++j ) {					
				for( index_type k=1; k < cfd.v.nZ-1 ; ++k ) {
					rt = cfd.iV2S(i,j,k,1);
					lt = cfd.iV2S(i,j,k,0);

					if( !cfd.isKinematicBoundary(lt) && !cfd.isKinematicBoundary(rt)) {
						if( phi(lt) > zero && phi(rt) < zero ) {
							theta = phi(rt) / (phi(rt) - phi(lt));
							inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

							v(i,j,k) += tension_factor * inv_theta * lerp(curvature3D(lt), curvature3D(rt), theta);
						}
						else if( phi(lt) < zero && phi(rt) > zero ) {
							theta = phi(lt) / (phi(lt) - phi(rt));
							inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

							v(i,j,k) -= tension_factor * inv_theta * lerp(curvature3D(rt), curvature3D(lt), theta);
						}
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
			index_type rt, lt;
			double theta, inv_theta;

			for( index_type j=1; j < cfd.w.nY-1 ; ++j ) {					
				for( index_type k=1; k < cfd.w.nZ-1 ; ++k ) {

					rt = cfd.iW2S(i,j,k,1);
					lt = cfd.iW2S(i,j,k,0);

					if( !cfd.isKinematicBoundary(lt) && !cfd.isKinematicBoundary(rt)) {
						if( phi(lt) > zero && phi(rt) < zero ) {
							theta = phi(rt) / (phi(rt) - phi(lt));
							inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

							w(i,j,k) += tension_factor * inv_theta * lerp(curvature3D(lt), curvature3D(rt), theta);
						}
						else if( phi(lt) < zero && phi(rt) > zero ) {
							theta = phi(lt) / (phi(lt) - phi(rt));
							inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

							w(i,j,k) -= tension_factor * inv_theta * lerp(curvature3D(rt), curvature3D(lt), theta);
						}
					}
				}			
			}
		}
	}

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

void CFDProjection_FreeSurface_Gibou::buildPoissonMatrix(Grid3d& u_part, Grid3d& v_part, Grid3d& w_part, Grid3d& s_part, const Grid3r& phi, const CFDsolver3D& cfd) {
	
	if(applySurfaceTension) {
		curvature3D.resize(cfd.isKinematicBoundary);
		tracking->getCurvature(curvature3D);

		log_info << "Apply Tension in Gibou Projection !!!!!" << endl;
	}

	if(u_part.is_empty())	u_part.resize(cfd.u);
	if(v_part.is_empty())	v_part.resize(cfd.v);
	if(w_part.is_empty())	w_part.resize(cfd.w);
	if(s_part.is_empty())	s_part.resize(cfd.isKinematicBoundary);

	const Grid3b&	isKinematicBoundary = cfd.isKinematicBoundary;
	const double	one = 1.0;	

	const real zero = real(0.0);

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
		double theta, inv_theta;

		for ( index_type j = 1; j < isKinematicBoundary.nY-1; ++j ) {
			for( index_type k = 1; k < isKinematicBoundary.nZ-1; ++k )	{
				

				if( phi(i,j,k) < zero && !isKinematicBoundary(i,j,k) )	{
					if( !isKinematicBoundary(i-1,j,k) && phi(i-1,j,k) >= zero )	{					
						theta = phi(i,j,k) / (phi(i,j,k) - phi(i-1,j,k));
						inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

						s_part(i,j,k) += inv_theta;
					}

					if( !isKinematicBoundary(i+1,j,k) )	{	
						if(phi(i+1,j,k) < zero) {
							u_part(i,j,k) += one;		
						}
						else {
							theta = phi(i,j,k) / (phi(i,j,k) - phi(i+1,j,k));
							inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

							s_part(i,j,k) += inv_theta;
						}	
					}

					if( !isKinematicBoundary(i,j-1,k) && phi(i,j-1,k) >= zero) {	
						theta = phi(i,j,k) / (phi(i,j,k) - phi(i,j-1,k));
						inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

						s_part(i,j,k) += inv_theta;
					}

					if( !isKinematicBoundary(i,j+1,k) )	{						
						if(phi(i,j+1,k) < zero) {
							v_part(i,j,k) += one;		
						}
						else {
							theta = phi(i,j,k) / (phi(i,j,k) - phi(i,j+1,k));
							inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

							s_part(i,j,k) += inv_theta;
						}	
					}

					if( !isKinematicBoundary(i,j,k-1) && phi(i,j,k-1) >= zero )	{					
						theta = phi(i,j,k) / (phi(i,j,k) - phi(i,j,k-1));
						inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

						s_part(i,j,k) += inv_theta;
					}

					if( !isKinematicBoundary(i,j,k+1) )	{						
						if(phi(i,j,k+1) < zero) {
							w_part(i,j,k) += one;		
						}
						else {
							theta = phi(i,j,k) / (phi(i,j,k) - phi(i,j,k+1));
							inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

							s_part(i,j,k) += inv_theta;
						}	
					}
				}
			}
		}
	}
}

void CFDProjection_FreeSurface_Gibou::addPressureGradient(CFDsolver2D& cfd, const Grid2d& pressure, const Grid2r& phi) {

	const real zero				= real(0.0);
	const real tension_factor	= sigma * dt * cfd.hInv;
	
	for( index_type i=1; i < cfd.u.nX-1 ; ++i ) {
		index_type rt, lt;
		double theta, inv_theta;

		for( index_type j=1; j < cfd.u.nY-1 ; ++j ) {					
			rt = cfd.iU2S(i,j,1);
			lt = cfd.iU2S(i,j,0);

			if( !cfd.isKinematicBoundary(lt) && !cfd.isKinematicBoundary(rt)) {
				if( phi(lt) < zero && phi(rt) < zero ) {
					cfd.u(i,j) += ( pressure(rt) - pressure(lt) ) ;
				}

				else if( phi(rt) < zero ) {			
					theta = phi(rt) / (phi(rt) - phi(lt));
					inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

					cfd.u(i,j) += ( pressure(rt) - zero ) * inv_theta ;
				}
				else if( phi(lt) < zero ) {
					theta = phi(lt) / (phi(lt) - phi(rt));
					inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

					cfd.u(i,j) += ( zero - pressure(lt) ) * inv_theta ;
				}
			}
		}
	}

	for( index_type i=1; i < cfd.v.nX-1 ; ++i ) {
		index_type rt, lt;
		double theta, inv_theta;

		for( index_type j=1; j < cfd.v.nY-1 ; ++j ) {					
			rt = cfd.iV2S(i,j,1);
			lt = cfd.iV2S(i,j,0);

			if( !cfd.isKinematicBoundary(lt) && !cfd.isKinematicBoundary(rt) ) {
				if( phi(lt) < zero && phi(rt) < zero ) {
					cfd.v(i,j) += ( pressure(rt) - pressure(lt) );
				}

				else if( phi(rt) < zero ) {
					theta = phi(rt) / (phi(rt) - phi(lt));
					inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

					cfd.v(i,j) += ( pressure(rt) - zero ) * inv_theta ;
				}
				else if( phi(lt) < zero ) {
					theta = phi(lt) / (phi(lt) - phi(rt));
					inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

					cfd.v(i,j) += ( zero - pressure(lt) ) * inv_theta ;
				}
			}
		}
	}
}

void CFDProjection_FreeSurface_Gibou::addPressureGradient(CFDsolver3D& cfd, const Grid3d& pressure, const Grid3r& phi) {	

	const real zero = real(0.0);

#ifdef USE_OPENMP
	int my_end = int(cfd.u.nX-1);
	int i;

#pragma omp parallel for
#else
	index_type my_end = cfd.u.nX-1;
	index_type i;
#endif
	for( i = 1; i < my_end; ++i ) {	
		index_type rt, lt;
		double theta, inv_theta;

		for( index_type j=1; j < cfd.u.nY-1 ; ++j ) {					
			for( index_type k=1; k < cfd.u.nZ-1 ; ++k ) {
				rt = cfd.iU2S(i,j,k,1);
				lt = cfd.iU2S(i,j,k,0);

				if( !cfd.isKinematicBoundary(lt) && !cfd.isKinematicBoundary(rt)) {
					if( phi(lt) < zero && phi(rt) < zero )
						cfd.u(i,j,k) += ( pressure(rt) - pressure(lt) ) ;

					else if( phi(rt) < zero ) {
						theta = phi(rt) / (phi(rt) - phi(lt));
						inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

						cfd.u(i,j,k) += ( pressure(rt) - zero ) * inv_theta;
					}
					else if( phi(lt) < zero ) {
						theta = phi(lt) / (phi(lt) - phi(rt));
						inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

						cfd.u(i,j,k) += ( zero - pressure(lt) ) * inv_theta;
					}
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
		index_type rt, lt;
		double theta, inv_theta;

		for( index_type j=1; j < cfd.v.nY-1 ; ++j ) {					
			for( index_type k=1; k < cfd.v.nZ-1 ; ++k ) {
				rt = cfd.iV2S(i,j,k,1);
				lt = cfd.iV2S(i,j,k,0);

				if( !cfd.isKinematicBoundary(lt) && !cfd.isKinematicBoundary(rt)) {
					if( phi(lt) < zero && phi(rt) < zero )
						cfd.v(i,j,k) += ( pressure(rt) - pressure(lt) );

					else if( phi(rt) < zero ) {
						theta = phi(rt) / (phi(rt) - phi(lt));
						inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

						cfd.v(i,j,k) += ( pressure(rt) - zero ) * inv_theta;
					}
					else if( phi(lt) < zero ) {
						theta = phi(lt) / (phi(lt) - phi(rt));
						inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

						cfd.v(i,j,k) += ( zero - pressure(lt) ) * inv_theta;
					}
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
		index_type rt, lt;
		double theta, inv_theta;

		for( index_type j=1; j < cfd.w.nY-1 ; ++j ) {					
			for( index_type k=1; k < cfd.w.nZ-1 ; ++k ) {

				rt = cfd.iW2S(i,j,k,1);
				lt = cfd.iW2S(i,j,k,0);

				if( !cfd.isKinematicBoundary(lt) && !cfd.isKinematicBoundary(rt)) {
					if( phi(lt) < zero && phi(rt) < zero )
						cfd.w(i,j,k) += ( pressure(rt) - pressure(lt) );

					else if( phi(rt) < zero ) {
						theta = phi(rt) / (phi(rt) - phi(lt));
						inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

						cfd.w(i,j,k) += ( pressure(rt) - zero ) * inv_theta;
					}
					else if( phi(lt) < zero ) {
						theta = phi(lt) / (phi(lt) - phi(rt));
						inv_theta = MAX(minInvTheta, MIN(maxInvTheta, (1.0f / theta)));

						cfd.w(i,j,k) += ( zero - pressure(lt) ) * inv_theta;
					}
				}
			}			
		}
	}
}