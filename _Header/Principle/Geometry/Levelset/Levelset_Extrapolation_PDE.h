#pragma once

#include "../../../Basis/Basis/Environments.h"

namespace geometryspace {

template<typename T>
void Levelset_Extrapolation_PDE<T>::extrapolate(Grid2<T,T>& val, typename MacGrid2D<T>::FIELD_TYPE _type, const GridBoundaryCondition<T,T>& bdry, const Grid2<T,T>& _phi) const {
	Grid2<T,T> phi; phi.resize(val);
	Levelset2D<T>::getPhi(phi, _phi, _type);
	extrapolate(val, bdry, phi);
}

template<typename T>
void Levelset_Extrapolation_PDE<T>::extrapolate(Grid3<T,T>& val, typename MacGrid3D<T>::FIELD_TYPE _type, const GridBoundaryCondition<T,T>& bdry, const Grid3<T,T>& _phi) const {
	Grid3<T,T> phi; phi.resize(val);
	Levelset3D<T>::getPhi(phi, _phi, _type);
	extrapolate(val, bdry, phi);
}

template<typename T>
void Levelset_Extrapolation_PDE<T>::extrapolate(Grid2<T,T>& val, const GridBoundaryCondition<T,T>& bdry, const Grid2<T,T>& phi) const {
	const T h = MAX(val.dx, val.dy);

	const T tau = cfl * h;
	const T eps = h, eps2 = eps*eps;
	const grid_index_type nIter = grid_index_type((narrowBandWidth+1) / cfl) + 1;

	Grid2<T,T> valTmp; valTmp.resize(val);

	//for( grid_index_type i = 1; i < val.nX-1; ++i ) {		
	//	for(grid_index_type j=1;j<val.nY-1;++j) {
	//		if(phi(i,j) > val.dx*3) {
	//			val(i,j) = 0;
	//			continue;
	//		}
	//	}
	//}

	bdry.enforce(val);
	for(grid_index_type iter = 0;iter < nIter; ++iter) {

#ifdef USE_OPENMP		
		int my_end = int(val.nX-1);
		int i;

#pragma omp parallel for
#else		
		grid_index_type my_end = val.nX-1;
		grid_index_type i;
#endif
		for( i = 1; i < my_end; ++i ) {		
			for(grid_index_type j=1;j<val.nY-1;++j) {

#define CLIP(x,xp,xm) if((x)>0) { xp = (x); xm = 0; } else { xp = 0; xm  = (x); }

				if(phi(i,j) <= 0) {
					valTmp(i,j) = val(i,j);
					continue;
				}

				const T _sign = phi(i,j) / sqrt(SQ(phi(i,j)) + eps2);

				const T a = val.dxInv * ( val(i,j) - val(i-1,j) );
				const T b = val.dxInv * ( val(i+1,j) - val(i,j) );
				const T c = val.dyInv * ( val(i,j) - val(i,j-1) );
				const T d = val.dyInv * ( val(i,j+1) - val(i,j) );				

				T phix = val.dxInv * ( phi(i+1,j) - phi(i-1,j) );
				T phiy = val.dyInv * ( phi(i,j+1) - phi(i,j-1) );
				T len = sqrt(phix*phix + phiy*phiy);

				Vector2<T> normal;
				if(IsAlmostZero(len)) {
					if(!IsAlmostSame(phi(i+1,j), phi(i,j)))	phix = val.dxInv * ( phi(i+1,j) - phi(i,j) );
					else									phix = val.dxInv * ( phi(i,j) - phi(i-1,j) );	

					if(!IsAlmostSame(phi(i,j+1), phi(i,j)))	phiy = val.dyInv * ( phi(i,j+1) - phi(i,j) );
					else									phiy = val.dyInv * ( phi(i,j) - phi(i,j-1) );	

					len = sqrt(phix*phix + phiy*phiy);

					if(IsAlmostZero(len))	normal = Vector2<T>(0,1);
					else					normal = Vector2<T>(phix / len, phiy / len);
				}
				else 
					normal = Vector2<T>(phix / len, phiy / len);					

				const T snx = _sign * normal[0];
				const T sny = _sign * normal[1];
				T snxp, snxm, snyp, snym;

				CLIP(snx, snxp, snxm);
				CLIP(sny, snyp, snym);

				valTmp(i,j) = val(i,j) - tau * (snxp * a + snxm * b + snyp * c + snym * d);				

#undef CLIP
			}
		}

		swap(val, valTmp);
		bdry.enforce(val);
	}	
}

template<typename T>
void Levelset_Extrapolation_PDE<T>::extrapolate(Grid3<T,T>& val, const GridBoundaryCondition<T,T>& bdry, const Grid3<T,T>& phi) const {
	const T h = MAX3(val.dx, val.dy, val.dz);

	const T tau = cfl * h;
	const T eps = h, eps2 = eps*eps;
	const grid_index_type nIter = grid_index_type((narrowBandWidth+1) / cfl) + 1;

	log_info << "Number of Iteration & Narrowband for extrapolation : " << nIter << " " << narrowBandWidth << std::endl; 

	Grid3<T,T> valTmp; valTmp.resize(val);

	bdry.enforce(val);
	for(grid_index_type iter = 0;iter < nIter; ++iter) {

#ifdef USE_OPENMP		
		int my_end = int(val.nX-1);
		int i;

#pragma omp parallel for
#else		
		grid_index_type my_end = val.nX-1;
		grid_index_type i;
#endif
		for( i = 1; i < my_end; ++i ) {		
			for(grid_index_type j=1;j<val.nY-1;++j) {
				for(grid_index_type k=1;k<val.nZ-1;++k) {

#define CLIP(x,xp,xm) if((x)>0) { xp = (x); xm = 0; } else { xp = 0; xm  = (x); }

					if(phi(i,j,k) <= 0) {
						valTmp(i,j,k) = val(i,j,k);
						continue;
					}

					T _sign = phi(i,j,k) / sqrt(SQ(phi(i,j,k)) + eps2);

					T a = val.dxInv * ( val(i,j,k) - val(i-1,j,k) );
					T b = val.dxInv * ( val(i+1,j,k) - val(i,j,k) );
					T c = val.dyInv * ( val(i,j,k) - val(i,j-1,k) );
					T d = val.dyInv * ( val(i,j+1,k) - val(i,j,k) );				
					T e = val.dzInv * ( val(i,j,k) - val(i,j,k-1) );
					T f = val.dzInv * ( val(i,j,k+1) - val(i,j,k) );				

					T phix = val.dxInv * ( phi(i+1,j,k) - phi(i-1,j,k) );
					T phiy = val.dyInv * ( phi(i,j+1,k) - phi(i,j-1,k) );
					T phiz = val.dzInv * ( phi(i,j,k+1) - phi(i,j,k-1) );

					T len = sqrt(phix*phix + phiy*phiy + phiz*phiz);

					Vector3<T> normal;
					if(IsAlmostZero(len)) {
						if(!IsAlmostSame(phi(i+1,j,k), phi(i,j,k)))	phix = val.dxInv * ( phi(i+1,j,k) - phi(i,j,k) );
						else										phix = val.dxInv * ( phi(i,j,k) - phi(i-1,j,k) );	

						if(!IsAlmostSame(phi(i,j+1,k), phi(i,j,k)))	phiy = val.dyInv * ( phi(i,j+1,k) - phi(i,j,k) );
						else										phiy = val.dyInv * ( phi(i,j,k) - phi(i,j-1,k) );	

						if(!IsAlmostSame(phi(i,j,k+1), phi(i,j,k)))	phiz = val.dzInv * ( phi(i,j,k+1) - phi(i,j,k) );
						else										phiz = val.dzInv * ( phi(i,j,k) - phi(i,j,k-1) );	

						len = sqrt(phix*phix + phiy*phiy + phiz*phiz);

						if(IsAlmostZero(len))	normal = Vector3<T>(0,1,0);
						else					normal = Vector3<T>(phix / len, phiy / len, phiz / len);
					}
					else 
						normal = Vector3<T>(phix / len, phiy / len, phiz / len);

					T snx = _sign * normal[0];
					T sny = _sign * normal[1];
					T snz = _sign * normal[2];
					T snxp, snxm, snyp, snym, snzp, snzm;

					CLIP(snx, snxp, snxm);
					CLIP(sny, snyp, snym);
					CLIP(snz, snzp, snzm);

					valTmp(i,j,k) = val(i,j,k) - tau * (snxp * a + snxm * b + snyp * c + snym * d + snzp * e + snzm * f);				

#undef CLIP

				}
			}
		}

		swap(val, valTmp);
		bdry.enforce(val);
	}	
}

};		// geometryspace