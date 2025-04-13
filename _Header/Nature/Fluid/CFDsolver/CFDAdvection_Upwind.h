#pragma once

#include "../../../Basis/Basis/Environments.h"

namespace fluidspace {

template<typename RK> template<typename DEV>
void CFDAdvection_Upwind<RK>::advect_ee(Grid2r& s, const Grid2r& s0, const CFDsolver2D& cfd, real dt, MacGrid2D<real>::FIELD_TYPE type, const GridBdry& b, real cfl) {
	Grid2r s_tmp(s0);

	const geometryspace::MacGrid2D_VelocityGetter<real> vel_getter(cfd.u, cfd.v, type);
	CFDTimeStepRestriction timeStepRestriction(cfl);

	const index_type sx_1 = s.nX-1;
	const index_type sy_1 = s.nY-1;

	real dtSum = 0;
	bool isEnd = false;

	while(!isEnd) {
		real _dt = dt - dtSum;
		real dt0 = _dt;
		
		if(cfl < 100)
			dt0 = timeStepRestriction.getAllowedTimeStep(cfd, MIN(s.dx, s.dy));

		if(dt0 < _dt)	_dt = dt0;
		else			isEnd = true;	

		dtSum += _dt;

#ifdef USE_OPENMP	
		int my_end = int(sx_1);
		int i;

#pragma omp parallel for
#else	
		index_type my_end = sx_1;
		index_type i;
#endif
		for( i = 1; i < my_end; ++i ) {		
			DEV upwind(s_tmp);
			arithmetic::Vector2<real> vel;
			real sx, sy;

			for(index_type j=1;j<sy_1;++j) {
				vel = vel_getter.getVel(i,j);

				sx = upwind.dx(i,j,vel[0] > 0);	
				sy = upwind.dy(i,j,vel[1] > 0);	

				s(i,j) = s_tmp(i,j) - _dt * ( vel[0] * sx + vel[1] * sy );
			}		
		}

		b.enforce(s);		
		swap(s, s_tmp);
	}	
	swap(s, s_tmp);
}	

template<typename RK> template<typename DEV>
void CFDAdvection_Upwind<RK>::advect_ee(Grid3r& s, const Grid3r& s0, const CFDsolver3D& cfd, real dt, MacGrid3D<real>::FIELD_TYPE type, const GridBdry& b, real cfl) {
	Grid3r s_tmp(s0);	

	const geometryspace::MacGrid3D_VelocityGetter<real> vel_getter(cfd.u, cfd.v, cfd.w, type);
	CFDTimeStepRestriction timeStepRestriction(cfl);

	const index_type sx_1 = s.nX-1;
	const index_type sy_1 = s.nY-1;
	const index_type sz_1 = s.nZ-1;

	real dtSum = 0;
	bool isEnd = false;

	while(!isEnd) {
		real _dt = dt - dtSum;
		real dt0 = _dt;

		if(cfl < 100)
			dt0 = timeStepRestriction.getAllowedTimeStep(cfd, MIN3(s.dx, s.dy, s.dz));

		if(dt0 < _dt)	_dt = dt0;
		else			isEnd = true;

		dtSum += _dt;

#ifdef USE_OPENMP	
		int my_end = int(sx_1);
		int i;

#pragma omp parallel for
#else	
		index_type my_end = sx_1;
		index_type i;
#endif
		for( i = 1; i < my_end; ++i ) {		
			DEV upwind(s_tmp);
			arithmetic::Vector3<real> vel;
			real sx, sy, sz;

			for(index_type j=1;j<sy_1;++j) {
				for(index_type k=1;k<sz_1;++k) {

					vel = vel_getter.getVel(i,j,k);

					sx = upwind.dx(i,j,k,vel[0] > 0);	
					sy = upwind.dy(i,j,k,vel[1] > 0);	
					sz = upwind.dz(i,j,k,vel[2] > 0);	

					s(i,j,k) = s_tmp(i,j,k) - _dt * ( vel[0] * sx + vel[1] * sy + vel[2] * sz );

				}			
			}		
		}

		b.enforce(s);
		swap(s, s_tmp);
	}	
	swap(s, s_tmp);
}

template<typename RK> template<typename DEV>
void CFDAdvection_Upwind<RK>::advect_ee(Grid2r& s, const Grid2r& s0, const VelocityField2D<real>& velfield, real dt, const GridBdry& b, real cfl) {
	
	Grid2r s_tmp(s0);

	CFDTimeStepRestriction timeStepRestriction(cfl);

	const index_type sx_1 = s.nX-1;
	const index_type sy_1 = s.nY-1;

	real dtSum = 0;
	bool isEnd = false;

	while(!isEnd) {
		real _dt = dt - dtSum;
		real dt0 = _dt;

		if(cfl < 100)
			dt0 = timeStepRestriction.getAllowedTimeStep(velfield, MIN(s.dx, s.dy));

		if(dt0 < _dt)	_dt = dt0;
		else			isEnd = true;		

		dtSum += _dt;

#ifdef USE_OPENMP	
		int my_end = int(sx_1);
		int i;

#pragma omp parallel for
#else	
		index_type my_end = sx_1;
		index_type i;
#endif

		for( i = 1; i < my_end; ++i ) {	
			DEV upwind(s_tmp);
			arithmetic::Vector2<real> vel;
			real sx, sy;

			for(index_type j=1;j<sy_1;++j) {
				vel = velfield.getVel(s.pos(i,j));

				sx = upwind.dx(i,j,vel[0] > 0);	
				sy = upwind.dy(i,j,vel[1] > 0);	

				s(i,j) = s_tmp(i,j) - _dt * ( vel[0] * sx + vel[1] * sy );
			}		
		}

		b.enforce(s);	
		swap(s, s_tmp);
	}	
	swap(s, s_tmp);
}

template<typename RK> template<typename DEV>
void CFDAdvection_Upwind<RK>::advect_ee(Grid3r& s, const Grid3r& s0, const VelocityField3D<real>& velfield, real dt, const GridBdry& b, real cfl) {

	Grid3r s_tmp(s0);

	CFDTimeStepRestriction timeStepRestriction(cfl);

	const index_type sx_1 = s.nX-1;
	const index_type sy_1 = s.nY-1;
	const index_type sz_1 = s.nZ-1;

	real dtSum = 0;
	bool isEnd = false;

	while(!isEnd) {
		real _dt = dt - dtSum;
		real dt0 = _dt;

		if(cfl < 100)
			dt0 = timeStepRestriction.getAllowedTimeStep(velfield, MIN3(s.dx, s.dy, s.dz));

		if(dt0 < _dt)	_dt = dt0;
		else			isEnd = true;	

		dtSum += _dt;

#ifdef USE_OPENMP	
		int my_end = int(sx_1);
		int i;

#pragma omp parallel for
#else	
		index_type my_end = sx_1;
		index_type i;
#endif

		for( i = 1; i < my_end; ++i ) {	
			DEV upwind(s_tmp);
			arithmetic::Vector3<real> vel;
			real sx, sy, sz;

			for(index_type j=1;j<sy_1;++j) {
				for(index_type k=1;k<sz_1;++k) {
					vel = velfield.getVel(s.pos(i,j,k));

					sx = upwind.dx(i,j,k,vel[0] > 0);	
					sy = upwind.dy(i,j,k,vel[1] > 0);	
					sz = upwind.dz(i,j,k,vel[2] > 0);	

					s(i,j,k) = s_tmp(i,j,k) - _dt * ( vel[0] * sx + vel[1] * sy + vel[2] * sz );
				}			
			}		
		}

		b.enforce(s);	
		swap(s, s_tmp);
	}	
	swap(s, s_tmp);
}

// 1st Upwind
template<typename RK> 
void CFDAdvection_1stUpwind<RK>::advect(Grid2r& s, const Grid2r& s0, const CFDsolver2D& cfd, real dt, MacGrid2D<real>::FIELD_TYPE type, const GridBdry& b) const {		
	RK()(s, s0, CFDAdvection_Upwind::Advect_EE_CFD<CFDsolver2D, Upwind_1st<Grid2r>>(cfd, dt, type, b, cfl));	
}

template<typename RK> 
void CFDAdvection_1stUpwind<RK>::advect(Grid3r& s, const Grid3r& s0, const CFDsolver3D& cfd, real dt, MacGrid3D<real>::FIELD_TYPE type, const GridBdry& b) const {
	RK()(s, s0, CFDAdvection_Upwind::Advect_EE_CFD<CFDsolver3D, Upwind_1st<Grid3r>>(cfd, dt, type, b, cfl));	
}

template<typename RK> 
void CFDAdvection_1stUpwind<RK>::advect(Grid2r& s, const Grid2r& s0, const VelocityField2D<real>& vel, real dt, const GridBdry& b) const {
	RK()(s, s0, CFDAdvection_Upwind::Advect_EE_VelField<VelocityField2D<real>, Upwind_1st<Grid2r>>(vel, dt, b, cfl));	
}

template<typename RK> 
void CFDAdvection_1stUpwind<RK>::advect(Grid3r& s, const Grid3r& s0, const VelocityField3D<real>& vel, real dt, const GridBdry& b) const {
	RK()(s, s0, CFDAdvection_Upwind::Advect_EE_VelField<VelocityField3D<real>, Upwind_1st<Grid3r>>(vel, dt, b, cfl));	
}	

// 2nd Upwind
template<typename RK> 
void CFDAdvection_2ndUpwind<RK>::advect(Grid2r& s, const Grid2r& s0, const CFDsolver2D& cfd, real dt, MacGrid2D<real>::FIELD_TYPE type, const GridBdry& b) const {		
	RK()(s, s0, CFDAdvection_Upwind::Advect_EE_CFD<CFDsolver2D, Upwind_2nd<Grid2r>>(cfd, dt, type, b, cfl));	
}

template<typename RK> 
void CFDAdvection_2ndUpwind<RK>::advect(Grid3r& s, const Grid3r& s0, const CFDsolver3D& cfd, real dt, MacGrid3D<real>::FIELD_TYPE type, const GridBdry& b) const {
	RK()(s, s0, CFDAdvection_Upwind::Advect_EE_CFD<CFDsolver3D, Upwind_2nd<Grid3r>>(cfd, dt, type, b, cfl));	
}

template<typename RK> 
void CFDAdvection_2ndUpwind<RK>::advect(Grid2r& s, const Grid2r& s0, const VelocityField2D<real>& vel, real dt, const GridBdry& b) const {
	RK()(s, s0, CFDAdvection_Upwind::Advect_EE_VelField<VelocityField2D<real>, Upwind_2nd<Grid2r>>(vel, dt, b, cfl));	
}

template<typename RK> 
void CFDAdvection_2ndUpwind<RK>::advect(Grid3r& s, const Grid3r& s0, const VelocityField3D<real>& vel, real dt, const GridBdry& b) const {
	RK()(s, s0, CFDAdvection_Upwind::Advect_EE_VelField<VelocityField3D<real>, Upwind_2nd<Grid3r>>(vel, dt, b, cfl));	
}	

// ENO
template<typename RK> 
void CFDAdvection_ENO<RK>::advect(Grid2r& s, const Grid2r& s0, const CFDsolver2D& cfd, real dt, MacGrid2D<real>::FIELD_TYPE type, const GridBdry& b) const {
	RK()(s, s0, CFDAdvection_Upwind::Advect_EE_CFD<CFDsolver2D, Upwind_ENO<Grid2r>>(cfd, dt, type, b, cfl));	
}

template<typename RK> 
void CFDAdvection_ENO<RK>::advect(Grid3r& s, const Grid3r& s0, const CFDsolver3D& cfd, real dt, MacGrid3D<real>::FIELD_TYPE type, const GridBdry& b) const {
	RK()(s, s0, CFDAdvection_Upwind::Advect_EE_CFD<CFDsolver3D, Upwind_ENO<Grid3r>>(cfd, dt, type, b, cfl));	
}

template<typename RK> 
void CFDAdvection_ENO<RK>::advect(Grid2r& s, const Grid2r& s0, const VelocityField2D<real>& vel, real dt, const GridBdry& b) const {
	RK()(s, s0, CFDAdvection_Upwind::Advect_EE_VelField<VelocityField2D<real>, Upwind_ENO<Grid2r>>(vel, dt, b, cfl));	
}

template<typename RK> 
void CFDAdvection_ENO<RK>::advect(Grid3r& s, const Grid3r& s0, const VelocityField3D<real>& vel, real dt, const GridBdry& b) const {
	RK()(s, s0, CFDAdvection_Upwind::Advect_EE_VelField<VelocityField3D<real>, Upwind_ENO<Grid3r>>(vel, dt, b, cfl));	
}

// WENO
template<typename RK> 
void CFDAdvection_WENO<RK>::advect(Grid2r& s, const Grid2r& s0, const CFDsolver2D& cfd, real dt, MacGrid2D<real>::FIELD_TYPE type, const GridBdry& b) const {
	RK()(s, s0, CFDAdvection_Upwind::Advect_EE_CFD<CFDsolver2D, Upwind_WENO<Grid2r>>(cfd, dt, type, b, cfl));	
}

template<typename RK> 
void CFDAdvection_WENO<RK>::advect(Grid3r& s, const Grid3r& s0, const CFDsolver3D& cfd, real dt, MacGrid3D<real>::FIELD_TYPE type, const GridBdry& b) const {
	RK()(s, s0, CFDAdvection_Upwind::Advect_EE_CFD<CFDsolver3D, Upwind_WENO<Grid3r>>(cfd, dt, type, b, cfl));	
}

template<typename RK> 
void CFDAdvection_WENO<RK>::advect(Grid2r& s, const Grid2r& s0, const VelocityField2D<real>& vel, real dt, const GridBdry& b) const {
	RK()(s, s0, CFDAdvection_Upwind::Advect_EE_VelField<VelocityField2D<real>, Upwind_WENO<Grid2r>>(vel, dt, b, cfl));	
}

template<typename RK> 
void CFDAdvection_WENO<RK>::advect(Grid3r& s, const Grid3r& s0, const VelocityField3D<real>& vel, real dt, const GridBdry& b) const {
	RK()(s, s0, CFDAdvection_Upwind::Advect_EE_VelField<VelocityField3D<real>, Upwind_WENO<Grid3r>>(vel, dt, b, cfl));	
}

};		// fluidspace