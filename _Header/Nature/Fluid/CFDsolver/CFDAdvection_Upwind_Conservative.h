#pragma once

namespace fluidspace {

template<class CFD> 
void CFDAdvection_Upwind_Conservative::advect_ee(Grid2r& s, const Grid2r& s0, const CFD& cfd, real dt, const GridBdry& b) const {
	Grid2r s_tmp(s0);

	CFDTimeStepRestriction timeStepRestriction(cfl);

	const index_type sx_1 = s.nX-1;
	const index_type sy_1 = s.nY-1;

	const real dx_half = s.dx * 0.5f, dxInv = (1.0 / s.dx);
	const real dy_half = s.dy * 0.5f, dyInv = (1.0 / s.dy);

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
			arithmetic::Vector2<real> pos;
			real u_0, u_1, v_0, v_1;

			for(index_type j=1;j<sy_1;++j) {
				pos = s.pos(i,j);

				u_0 = cfd.getU(pos - Vector2<real>(dx_half, 0));
				u_1 = cfd.getU(pos + Vector2<real>(dx_half, 0));
				v_0 = cfd.getV(pos - Vector2<real>(0, dy_half));
				v_1 = cfd.getV(pos + Vector2<real>(0, dy_half));

				s(i,j) = s_tmp(i,j) + _dt *	(( u_0 * ((u_0 > 0) ? s_tmp(i-1,j) : s_tmp(i  ,j))
											 - u_1 * ((u_1 > 0) ? s_tmp(i  ,j) : s_tmp(i+1,j)) ) * dxInv
										   + ( v_0 * ((v_0 > 0) ? s_tmp(i,j-1) : s_tmp(i,j  ))
											 - v_1 * ((v_1 > 0) ? s_tmp(i,j  ) : s_tmp(i,j+1)) ) * dyInv );
			}	
		}

		b.enforce(s);		
		swap(s, s_tmp);
	}	
	swap(s, s_tmp);
}

template<class CFD> 
void CFDAdvection_Upwind_Conservative::advect_ee(Grid3r& s, const Grid3r& s0, const CFD& cfd, real dt, const GridBdry& b) const {
	Grid3r s_tmp(s0);

	CFDTimeStepRestriction timeStepRestriction(cfl);

	const index_type sx_1 = s.nX-1;
	const index_type sy_1 = s.nY-1;
	const index_type sz_1 = s.nZ-1;

	const real dx_half = s.dx * 0.5f, dxInv = 1.0 / s.dx;
	const real dy_half = s.dy * 0.5f, dyInv = 1.0 / s.dy;
	const real dz_half = s.dz * 0.5f, dzInv = 1.0 / s.dz;

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
			arithmetic::Vector3<real> pos;
			real u_0, u_1, v_0, v_1, w_0, w_1;

			for(index_type j=1;j<sy_1;++j) {
				for(index_type k=1;k<sz_1;++k) {
					pos = s.pos(i,j,k);

					u_0 = cfd.getU(pos - Vector3<real>(dx_half, 0, 0));
					u_1 = cfd.getU(pos + Vector3<real>(dx_half, 0, 0));
					v_0 = cfd.getV(pos - Vector3<real>(0, dy_half, 0));
					v_1 = cfd.getV(pos + Vector3<real>(0, dy_half, 0));
					w_0 = cfd.getW(pos - Vector3<real>(0, 0, dz_half));
					w_1 = cfd.getW(pos + Vector3<real>(0, 0, dz_half));

					s(i,j,k) = s_tmp(i,j,k) + _dt *	( ( u_0 * ((u_0 > 0) ? s_tmp(i-1,j,k) : s_tmp(i  ,j,k))
													  - u_1 * ((u_1 > 0) ? s_tmp(i  ,j,k) : s_tmp(i+1,j,k)) ) * dxInv
													+ ( v_0 * ((v_0 > 0) ? s_tmp(i,j-1,k) : s_tmp(i,j  ,k))
													  - v_1 * ((v_1 > 0) ? s_tmp(i,j  ,k) : s_tmp(i,j+1,k)) ) * dyInv
													+ ( w_0 * ((w_0 > 0) ? s_tmp(i,j,k-1) : s_tmp(i,j,k  ))
													  - w_1 * ((w_1 > 0) ? s_tmp(i,j,k  ) : s_tmp(i,j,k+1)) ) * dzInv );
				}				
			}	
		}

		b.enforce(s);		
		swap(s, s_tmp);
	}	
	swap(s, s_tmp);
}

};		// fluidspace