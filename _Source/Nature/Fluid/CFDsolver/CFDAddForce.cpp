#include <stdafx.h>
#include <Nature/Fluid/CFDsolver/CFDAddForce.h>

using namespace std;
using namespace fluidspace;

void CFDAddForce::addForce(CFDsolver2D& cfd, const Grid2<Vector2<real>,real>& field, real dt) {
	const real half_dt = real(dt * 0.5);

	for(index_type i=1;i<cfd.u.nX-1;++i) 
		for(index_type j=1;j<cfd.u.nY-1;++j) 
			cfd.u(i,j) += ( field(i,j)[0] + field(i+1,j)[0] ) * half_dt;

	for(index_type i=1;i<cfd.v.nX-1;++i) 
		for(index_type j=1;j<cfd.v.nY-1;++j) 
			cfd.v(i,j) += ( field(i,j)[1] + field(i,j+1)[1] ) * half_dt;
}

void CFDAddForce::addForce(CFDsolver3D& cfd, const Grid3<Vector3<real>,real>& field, real dt) {
	const real half_dt = real(dt * 0.5);

#ifdef USE_OPENMP
	int my_end = int(cfd.u.nX-1);
	int i;

#pragma omp parallel for
#else
	index_type my_end = cfd.u.nX-1;
	index_type i;
#endif
	for(i=1;i<my_end;++i) 	
		for(index_type j=1;j<cfd.u.nY-1;++j) 
			for(index_type k=1;k<cfd.u.nZ-1;++k) 
				cfd.u(i,j,k) += ( field(i,j,k)[0] + field(i+1,j,k)[0] ) * half_dt;

#ifdef USE_OPENMP
	my_end = int(cfd.v.nX-1);
#pragma omp parallel for
#else
	my_end = cfd.v.nX-1;
#endif
	for(i=1;i<my_end;++i) 	
		for(index_type j=1;j<cfd.v.nY-1;++j) 
			for(index_type k=1;k<cfd.v.nZ-1;++k) 
				cfd.v(i,j,k) += ( field(i,j,k)[1] + field(i,j+1,k)[1] ) * half_dt;

#ifdef USE_OPENMP
	my_end = int(cfd.w.nX-1);
#pragma omp parallel for
#else
	my_end = cfd.w.nX-1;
#endif
	for(i=1;i<my_end;++i) 	
		for(index_type j=1;j<cfd.w.nY-1;++j) 
			for(index_type k=1;k<cfd.w.nZ-1;++k) 
				cfd.w(i,j,k) += ( field(i,j,k)[2] + field(i,j,k+1)[2] ) * half_dt;
}

void CFDAddForce_Gravity::addForce(CFDsolver2D& cfd, real dt) {
	const Grid2b&	isbdry = cfd.isKinematicBoundary;
	Grid2r&			v = cfd.v;

	for(index_type i=1;i<v.nX-1;++i) {
		for(index_type j=1;j<v.nY-1;++j) {
			if(!isbdry(i,j) && !isbdry(i,j+1))
				v(i,j) += dt * g;
		}
	}
}

void CFDAddForce_Gravity::addForce(CFDsolver3D& cfd, real dt) {

	const Grid3b&	isbdry = cfd.isKinematicBoundary;
	Grid3r&			v = cfd.v;

#ifdef USE_OPENMP
	int my_end = int(v.nX-1);
	int i;

#pragma omp parallel for
#else
	index_type my_end = v.nX-1;
	index_type i;
#endif
	for(i=1;i<my_end;++i) 
		for(index_type j=1;j<v.nY-1;++j) 
			for(index_type k=1;k<v.nZ-1;++k) 
				if(!isbdry(i,j,k) && !isbdry(i,j+1,k))
					v(i,j,k) += dt * g;
}

void CFDAddForce_DensityForce::addForce(CFDsolver2D& cfd, real dt) {
	if(density == NULL || density->getDimension() != 2)
		basis::throwError("Error in CFDAddForce_DensityForce::addForce(CFDsolver2D& cfd, real dt)");

	const Grid2b&	isbdry = cfd.isKinematicBoundary;
	Grid2r&			v = cfd.v;

	const real		factor = dt * (g + buoyancy);
	const Grid2r&	field = static_cast<const DensityTracking2D*>(density)->density;

	for( index_type i=1; i < v.nX-1; ++i ) 
		for( index_type j=1; j < v.nY-1; ++j ) 
			if(!isbdry(i,j) && !isbdry(i,j+1))
				v(i,j) += (factor * field(v.pos(i,j)));
}

void CFDAddForce_DensityForce::addForce(CFDsolver3D& cfd, real dt) {
	if(density == NULL || density->getDimension() != 3)
		basis::throwError("Error in CFDAddForce_DensityForce::addForce(CFDsolver3D& cfd, real dt)");

	const Grid3b&	isbdry = cfd.isKinematicBoundary;
	Grid3r&			v = cfd.v;

	const real		factor = dt * (g + buoyancy);
	const Grid3r&	field = static_cast<const DensityTracking3D*>(density)->getDensity();

#ifdef USE_OPENMP
	int my_end = int(v.nX-1);
	int i;

#pragma omp parallel for
#else
	index_type my_end = v.nX-1;
	index_type i;
#endif
	for( i = 1 ; i < my_end; ++i ) {
		for( index_type j=1; j < v.nY-1; ++j ) {
			for( index_type k=1; k < v.nZ-1; ++k ) {
				if(!isbdry(i,j,k) && !isbdry(i,j+1,k))
					v(i,j,k) += (factor * field(v.pos(i,j,k)));
			}}}
}

void CFDAddForce_SurfaceTension::addForce(CFDsolver2D& cfd, real dt) {
	CFDAddForce_Gravity::addForce(cfd, dt);

	if(tracking == NULL)
		return;

	Grid2r phi;	
	Grid2<Vector2<real>,real> forces; 

	phi.resize(cfd.isKinematicBoundary);
	forces.resize(cfd.isKinematicBoundary);

	tracking->getPhi(phi);	
	forces.zero();
	
	for(index_type i=1;i<forces.nX-1;++i) {
		real _delta, _curvature, _phi, _density;
		Vector2<real> _normal;

		for(index_type j=1;j<forces.nY-1;++j) {

			_phi = phi(i,j) * phi.dxInv;
			_delta = MacGridLibrary<real>::delta(_phi);

			if(_delta == 0)
				continue;

			_curvature = MacGridLibrary<real>::getCurvature(phi, i, j);
			_normal = MacGridLibrary<real>::getNormal(phi, i, j);
			_density = MacGridLibrary<real>::heaviside(_phi, rho_pos, rho_neg);

			forces(i,j) -= ((_delta * sigma * _curvature / _density) * _normal);

		}
	}	

	CFDAddForce::addForce(cfd, forces, dt);
}

void CFDAddForce_SurfaceTension::addForce(CFDsolver3D& cfd, real dt) {
	CFDAddForce_Gravity::addForce(cfd, dt);

	if(tracking == NULL)
		return;

	Grid3r phi;	
	Grid3<Vector3<real>,real> forces; 

	phi.resize(cfd.isKinematicBoundary);
	forces.resize(cfd.isKinematicBoundary);

	tracking->getPhi(phi);	
	forces.zero();

#ifdef USE_OPENMP
	int my_end = int(forces.nX-1);
	int i;

#pragma omp parallel for
#else
	index_type my_end = forces.nX-1;
	index_type i;
#endif
	for(i=1;i<my_end;++i) {
		real _delta, _curvature, _phi, _density;
		Vector3<real> _normal;

		for(index_type j=1;j<forces.nY-1;++j) {
			for(index_type k=1;k<forces.nZ-1;++k) {

				_phi = phi(i,j,k) * phi.dxInv;
				_delta = MacGridLibrary<real>::delta(_phi);

				if(_delta == 0)
					continue;

				_curvature = MacGridLibrary<real>::getCurvature(phi, i, j, k);
				_normal = MacGridLibrary<real>::getNormal(phi, i, j, k);
				_density = MacGridLibrary<real>::heaviside(_phi, rho_pos, rho_neg);

				forces(i,j,k) -= ((_delta * sigma * _curvature / _density) * _normal);

			}
		}
	}

	CFDAddForce::addForce(cfd, forces, dt);
}

//void CFDAddForce_Liquid_Bouyancy::addForce(CFDsolver2D& cfd, real dt) {
//	Grid2r phi;
//	tracking->getPhi(phi);
//
//#ifdef USE_OPENMP
//	int my_end = int(cfd.v.nX);
//	int i;
//
//#pragma omp parallel for
//#else
//	index_type my_end = cfd.v.nX;
//	index_type i;
//#endif
//	for( i = 0; i < my_end; ++i ) {
//		for( index_type j = 0; j < cfd.v.nY; ++j ) {
//			real _phi = phi.get(cfd.v.pos(i,j));
//			cfd.v(i,j) += dt * g * MacGridLibrary<real>::heaviside(_phi * cfd.v.dx, rho_pos, rho_neg);	
//		}		
//	}
//}
//
//void CFDAddForce_Liquid_Bouyancy::addForce(CFDsolver3D&, real dt) {
//	basis::throwError("Not Implemented _ CFDAddForce_Liquid_Bouyancy::addForce(CFDsolver3D&, real dt)");
//}