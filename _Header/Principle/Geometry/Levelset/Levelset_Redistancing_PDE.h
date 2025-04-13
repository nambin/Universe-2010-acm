#pragma once

namespace geometryspace {

template<typename T, typename RK>
const T Levelset_Redistancing_PDE<T,RK>::sign_sussman94(const T phi, const T h_sqr) {
	const T _phi = (phi == 0) ? -1.0e-6 : phi;
	return T(_phi / sqrt(_phi*_phi + h_sqr));
}

template<typename T, typename RK>
const T Levelset_Redistancing_PDE<T,RK>::sign_peng99(const Grid2<T,T>& phi, const index_type i, const index_type j) {
	const T h = phi.dx, h2 = h*h, eps = 1.0e-6;
	const T _phi = (phi(i,j) == 0) ? -eps : phi(i,j);

	T denom2 = SQ(_phi) + (h2 * T(0.25)) * ( SQ( phi.dxInv * (phi(i+1,j) - phi(i-1,j)) ) + SQ( phi.dyInv * (phi(i,j+1) - phi(i,j-1)) ) );

	return T(_phi / sqrt(denom2));
}

template<typename T, typename RK>
const T Levelset_Redistancing_PDE<T,RK>::sign_peng99(const Grid3<T,T>& phi, const index_type i, const index_type j, const index_type k) {
	const T h = phi.dx, h2 = h*h, eps = 1.0e-6;
	const T _phi = (phi(i,j,k) == 0) ? -eps : phi(i,j,k);

	T denom2 = SQ(_phi) + (h2 * T(0.25)) * ( SQ( phi.dxInv * (phi(i+1,j,k) - phi(i-1,j,k)) ) + SQ( phi.dyInv * (phi(i,j+1,k) - phi(i,j-1,k)) ) + SQ( phi.dzInv * (phi(i,j,k+1) - phi(i,j,k-1)) ) );

	return T(_phi / sqrt(denom2));
}

template<typename T, typename RK> template<typename SWEEP> 
void Levelset_Redistancing_PDE<T,RK>::redistance_pde(Grid2<T,T>& phi, const Grid2<T,T>& phi0, const SWEEP& one_sweep) const {
	const T tau = cfl * MIN(phi.dx, phi.dy);
	grid_index_type nIter = grid_index_type((narrowBandWidth+1) / cfl) + 1;

	nIter = MAX((grid_index_type)min_iter, nIter);

	Grid2<T,T> phiTmp; phiTmp.resize(phi);
	if(&phi != &phi0)
		phi = phi0;

	one_sweep.set_tau(tau);
	for(grid_index_type iter = 0;iter < nIter; ++iter) {
		RK()(phiTmp, phi, one_sweep);
		swap(phi, phiTmp);
		phi.setNeumannBoundary();	
	}	
}

template<typename T, typename RK> template<typename SWEEP> 
void Levelset_Redistancing_PDE<T,RK>::redistance_pde(Grid3<T,T>& phi, const Grid3<T,T>& phi0, const SWEEP& one_sweep) const {
	const T tau = cfl * MIN3(phi.dx, phi.dy, phi.dz);
	grid_index_type nIter = grid_index_type((narrowBandWidth+1) / cfl) + 1;

	nIter = MAX((grid_index_type)min_iter, nIter);

	log_info << "Number of Iteration & Narrowband for redistancing levelset : " << nIter << " " << narrowBandWidth << std::endl; 

	Grid3<T,T> phiTmp; phiTmp.resize(phi);
	if(&phi != &phi0)
		phi = phi0;

	one_sweep.set_tau(tau);
	for(grid_index_type iter = 0;iter < nIter; ++iter) {
		RK()(phiTmp, phi, one_sweep);
		swap(phi, phiTmp);
		phi.setNeumannBoundary();	
	}	
}

template<typename T, typename RK> template<typename DEV> 
void Levelset_Redistancing_PDE<T,RK>::redistance_pde_one_sweep(Grid2<T,T>& phiTmp, const Grid2<T,T>& phi, const T tau) const {

	const T h = phi.dx, h2 = h*h;

#ifdef USE_OPENMP		
	int my_end = int(phi.nX-1);
	int i;

#pragma omp parallel for
#else		
	grid_index_type my_end = phi.nX-1;
	grid_index_type i;
#endif
	for( i = 1; i < my_end; ++i ) {	

		T a,b,c,d,sign;		
		T ap,bp,cp,dp;
		T am,bm,cm,dm;
		DEV upwind(phi);

		for(grid_index_type j=1;j<phi.nY-1;++j) {

#define CLIP(x,xp,xm) if((x)>0) { xp = (x); xm = 0; } else { xp = 0; xm  = (x); }

			upwind.dx(a,b,i,j);
			upwind.dy(c,d,i,j);
			
			CLIP(a,ap,am);	CLIP(b,bp,bm);
			CLIP(c,cp,cm);	CLIP(d,dp,dm);

			if(use_sign_peng)	sign = sign_peng99(phi, i, j);
			else				sign = sign_sussman94(phi(i,j), h2);

			if(sign > 0) {

				//	phiTmp(i,j) = phi(i,j) 
				//				- tau * sign * (sqrt(MAX(ap*ap, bm*bm) + MIN(cp*cp,dm*dm)) - 1);

				phiTmp(i,j) = phi(i,j) 
							- tau * sign * (sqrt(ap*ap + bm*bm + cp*cp + dm*dm) - 1);
			}
			else {

				//	phiTmp(i,j) = phi(i,j) 
				//				- tau * sign * (sqrt(MAX(am*am, bp*bp) + MIN(cm*cm, dp*dp)) - 1);

				phiTmp(i,j) = phi(i,j) 
							- tau * sign * (sqrt(am*am + bp*bp + cm*cm + dp*dp) - 1);

			}

#undef CLIP
		}
	}	
}

template<typename T, typename RK> template<typename DEV> 
void Levelset_Redistancing_PDE<T,RK>::redistance_pde_one_sweep(Grid3<T,T>& phiTmp, const Grid3<T,T>& phi, const T tau) const {
	const T h = phi.dx, h2 = h*h;

	const T eps = 1.0e-5 * h;

#ifdef USE_OPENMP		
	int my_end = int(phi.nX-1);
	int i;

#pragma omp parallel for
#else		
	grid_index_type my_end = phi.nX-1;
	grid_index_type i;
#endif
	for( i = 1; i < my_end; ++i ) {	

		T a,b,c,d,e,f,sign;		
		T ap,bp,cp,dp,ep,fp;
		T am,bm,cm,dm,em,fm;

		DEV upwind(phi);

		for(grid_index_type j=1;j<phi.nY-1;++j) {
			for(grid_index_type k=1;k<phi.nZ-1;++k) {

#define CLIP(x,xp,xm) if((x)>0) { xp = (x); xm = 0; } else { xp = 0; xm  = (x); }

				upwind.dx(a,b,i,j,k);
				upwind.dy(c,d,i,j,k);
				upwind.dz(e,f,i,j,k);

				CLIP(a,ap,am);	CLIP(b,bp,bm);
				CLIP(c,cp,cm);	CLIP(d,dp,dm);
				CLIP(e,ep,em);	CLIP(f,fp,fm);

				if(use_sign_peng)	sign = sign_peng99(phi, i, j, k);
				else				sign = sign_sussman94(phi(i,j,k), h2);

				if(sign > 0) {
					phiTmp(i,j,k) = phi(i,j,k) 
									- tau * sign * (sqrt(ap*ap + bm*bm + cp*cp + dm*dm + ep*ep + fm*fm) - 1);

					if(phiTmp(i,j,k) < 0) phiTmp(i,j,k) = + eps;
				}
				else {
					phiTmp(i,j,k) = phi(i,j,k) 
									- tau * sign * (sqrt(am*am + bp*bp + cm*cm + dp*dp + em*em + fp*fp) - 1);

					if(phiTmp(i,j,k) > 0) phiTmp(i,j,k) = - eps;
				}
#undef CLIP
			}
		}
	}	
}

};			// geometryspace