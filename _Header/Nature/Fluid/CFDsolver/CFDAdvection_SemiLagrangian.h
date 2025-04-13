#pragma once

#include "../../../Basis/Basis/Environments.h"
#include "../../../Principle/Geometry/MacGrid/MacGridFunctor.h"
#include "../CFDsolver/CFDAdvection.h"

namespace fluidspace {

template<typename T, typename TRACER> 
class CFD2D_advect_SemiLagrangian {
	Grid2<T,T>						& s;
	const Grid2<T,T>				& s0;

	const CFDsolver2D				& cfd;
	MacGrid2D<real>::FIELD_TYPE		type;

	T								dt; 

public : 
	typedef typename Grid2<T,T>::index_type index_type;

	CFD2D_advect_SemiLagrangian(Grid2<T,T> &_s, const Grid2<T,T> &_s0, const CFDsolver2D& _cfd, MacGrid2D<real>::FIELD_TYPE _type, T _dt)
	: s(_s), s0(_s0), cfd(_cfd), type(_type), dt(_dt) {}

	void operator()(const basis::blocked_range2d<index_type>& r) const {				
		Grid2b isbdry; isbdry.resize(s);
		cfd.interp_KinematicBoundary(isbdry, type);

		const AABB2D<T> aabb = cfd.getAABB();		

#ifdef USE_OPENMP
		int my_begin = int(r.rows().begin());
		int my_end = int(r.rows().end());
		int i;

#pragma omp parallel for
#else
		index_type my_begin = r.rows().begin();
		index_type my_end = r.rows().end();
		index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {			
			Vector2<T> posNew, path;

			for( index_type j=r.cols().begin(); j < r.cols().end(); ++j ) {
				if(isbdry(i,j)) {
					continue;					
				}

				posNew = TRACER::backtrack(s.pos(i,j), cfd, aabb, dt);
				cfd.constrainedBdry(posNew);
				
				s(i,j) = s0.get(posNew);
			}
		}
	}
};

template<typename T, typename TRACER> 
class CFD3D_advect_SemiLagrangian {
	Grid3<T,T>						& s;
	const Grid3<T,T>				& s0;

	const CFDsolver3D				& cfd;
	MacGrid3D<real>::FIELD_TYPE		type;

	T								dt; 

public : 
	typedef typename Grid3<T,T>::index_type index_type;

	CFD3D_advect_SemiLagrangian(Grid3<T,T> &_s, const Grid3<T,T> &_s0, const CFDsolver3D& _cfd, MacGrid3D<real>::FIELD_TYPE _type, T _dt)
	: s(_s), s0(_s0), cfd(_cfd), type(_type), dt(_dt) {}

	void operator()(const basis::blocked_range3d<index_type>& r) const {	

		Grid3b isbdry; isbdry.resize(s);
		cfd.interp_KinematicBoundary(isbdry, type);

		const AABB<T> aabb = cfd.getAABB();

#ifdef USE_OPENMP
		int my_begin = int(r.pages().begin());
		int my_end = int(r.pages().end());
		int i;

#pragma omp parallel for
#else
		index_type my_begin = r.pages().begin();
		index_type my_end = r.pages().end();
		index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {		
			Vector3<T> posNew;

			for( index_type j=r.rows().begin(); j < r.rows().end(); ++j ) {
				for( index_type k=r.cols().begin(); k < r.cols().end(); ++k ) {
					if(isbdry(i,j,k))
						continue;

					posNew = TRACER::backtrack(s.pos(i,j,k), cfd, aabb, dt);
					cfd.constrainedBdry(posNew);

					s(i,j,k) = s0.get(posNew);
				}
			}
		}
	}
};

template<typename TRACER>
void CFDAdvection_SemiLagrangian<TRACER>::advect_sl(Grid2r& s, const Grid2r& s0, const CFDsolver2D& cfd, real dt, MacGrid2D<real>::FIELD_TYPE type, const GridBdry& b) {
	CFD2D_advect_SemiLagrangian<real,TRACER> advector(s, s0, cfd, type, dt);
	advector(basis::blocked_range2d<index_type>(1, s.nX-1, 1, s.nY-1));

	b.enforce(s);
}

template<typename TRACER>
void CFDAdvection_SemiLagrangian<TRACER>::advect_sl(Grid3r& s, const Grid3r& s0, const CFDsolver3D& cfd, real dt, MacGrid3D<real>::FIELD_TYPE type, const GridBdry& b) {
	CFD3D_advect_SemiLagrangian<real,TRACER> advector(s, s0, cfd, type, dt);
	advector(basis::blocked_range3d<index_type>(1, s.nX-1, 1, s.nY-1, 1, s.nZ-1));

	b.enforce(s);
}

template<typename TRACER>
void CFDAdvection_SemiLagrangian<TRACER>::advect_sl(Grid2r& s, const Grid2r& s0, const VelocityField2D<real>& vel, real dt, const GridBdry& b) {
	VelocityField2D_advect_SemiLagrangian<real,VelocityField2D<real>,TRACER> advector(s, s0, vel, dt);
	advector(basis::blocked_range2d<index_type>(1, s.nX-1, 1, s.nY-1));

	b.enforce(s);
}	

template<typename TRACER>
void CFDAdvection_SemiLagrangian<TRACER>::advect_sl(Grid3r& s, const Grid3r& s0, const VelocityField3D<real>& vel, real dt, const GridBdry& b) {
	VelocityField3D_advect_SemiLagrangian<real,VelocityField3D<real>,TRACER> advector(s, s0, vel, dt);
	advector(basis::blocked_range3d<index_type>(1, s.nX-1, 1, s.nY-1, 1, s.nZ-1));

	b.enforce(s);
}

template<typename TRACER>
void CFDAdvection_SemiLagrangian<TRACER>::limiting_sl(Grid2r& s, const Grid2r& s0, const CFDsolver2D& cfd, real dt, MacGrid2D<real>::FIELD_TYPE type) {
	
	const MacGrid2D_VelocityGetter<real> vel(cfd.u, cfd.v, type);
	
#ifdef USE_OPENMP
	int my_end = int(s.nX-1);
	int i;

#pragma omp parallel for
#else
	index_type my_end = s.nX-1;
	index_type i;
#endif
	for( i = 1; i < my_end; ++i ) {		
		Vector2<real> posNew;
		real sMax, sMin;

		for ( index_type j = 1; j < s.nY-1; ++j ) {
			posNew = TRACER::backtrack(s.pos(i,j), vel, cfd, dt);
			cfd.constrainedBdry(posNew);

			s0.getMaxAndMin(sMax, sMin, posNew);			

			if ( s(i,j) > sMax )		s(i,j) = sMax;
			else if ( s(i,j) < sMin )	s(i,j) = sMin;
		}
	}	
}

template<typename TRACER>
void CFDAdvection_SemiLagrangian<TRACER>::limiting_sl(Grid2r& s, const Grid2r& s0, const VelocityField2D<real>& vel, real dt) {

	const AABB2D<real> aabb = s.getAABB();

#ifdef USE_OPENMP
	int my_end = int(s.nX-1);
	int i;

#pragma omp parallel for
#else
	index_type my_end = s.nX-1;
	index_type i;
#endif
	for( i = 1; i < my_end; ++i ) {			
		Vector2<real> posNew;
		real sMax, sMin;

		for ( index_type j = 1; j < s.nY-1; ++j ) {
			posNew = TRACER::backtrack(s.pos(i,j), vel, aabb, dt);
			aabb.constrainedBdry(posNew);

			s0.getMaxAndMin(sMax, sMin, posNew);

			if ( s(i,j) > sMax )		s(i,j) = sMax;
			else if ( s(i,j) < sMin )	s(i,j) = sMin;
		}
	}	
}

template<typename TRACER>
void CFDAdvection_SemiLagrangian<TRACER>::limiting_sl(Grid3r& s, const Grid3r& s0, const CFDsolver3D& cfd, real dt, MacGrid3D<real>::FIELD_TYPE type) {
	const MacGrid3D_VelocityGetter<real> vel(cfd.u, cfd.v, cfd.w, type);

#ifdef USE_OPENMP
	int my_end = int(s.nX-1);
	int i;

#pragma omp parallel for
#else
	index_type my_end = s.nX-1;
	index_type i;
#endif
	for( i = 1; i < my_end; ++i ) {		
		Vector3<real> posNew;
		real sMax, sMin;

		for ( index_type j = 1; j < s.nY-1; ++j ) {
			for ( index_type k = 1; k < s.nZ-1; ++k ) {

				posNew = TRACER::backtrack(s.pos(i,j,k), vel, cfd, dt);
				cfd.constrainedBdry(posNew);

				s0.getMaxAndMin(sMax, sMin, posNew);			

				if ( s(i,j,k) > sMax )		s(i,j,k) = sMax;
				else if ( s(i,j,k) < sMin )	s(i,j,k) = sMin;

			}			
		}
	}	
}

template<typename TRACER>
void CFDAdvection_SemiLagrangian<TRACER>::limiting_sl(Grid3r& s, const Grid3r& s0, const VelocityField3D<real>& vel, real dt) {
	const AABB<real> aabb = s.getAABB();

#ifdef USE_OPENMP
	int my_end = int(s.nX-1);
	int i;

#pragma omp parallel for
#else
	index_type my_end = s.nX-1;
	index_type i;
#endif
	for( i = 1; i < my_end; ++i ) {			
		Vector3<real> posNew;
		real sMax, sMin;

		for ( index_type j = 1; j < s.nY-1; ++j ) {
			for ( index_type k = 1; k < s.nZ-1; ++k ) {

				posNew = TRACER::backtrack(s.pos(i,j,k), vel, aabb, dt);
				aabb.constrainedBdry(posNew);

				s0.getMaxAndMin(sMax, sMin, posNew);			

				if ( s(i,j,k) > sMax )		s(i,j,k) = sMax;
				else if ( s(i,j,k) < sMin )	s(i,j,k) = sMin;

			}			
		}
	}	
}

template<typename TRACER>
void CFDAdvection_SemiLagrangian<TRACER>::revert_sl(Grid2r& s, const Grid2r& s0, const Grid2r& sTmp, const CFDsolver2D& cfd, real dt, MacGrid2D<real>::FIELD_TYPE type) {
	Vector2<real> posNew;
	real sMax, sMin;

	MacGrid2D_VelocityGetter<real> vel(cfd.u, cfd.v, type);

	// Limiting
	for ( index_type i = 1; i < s.nX-1; ++i ) {
		for ( index_type j = 1; j < s.nY-1; ++j ) {
			posNew = TRACER::backtrack(s.pos(i,j), vel, cfd, dt);
			cfd.constrainedBdry(posNew);

			s0.getMaxAndMin(sMax, sMin, posNew);

			if ( sMin < sTmp(i,j) && sTmp(i,j) < sMax ) s(i,j) = sTmp(i,j);			
		}
	}
}

template<typename TRACER>
void CFDAdvection_SemiLagrangian<TRACER>::revert_sl(Grid2r& s, const Grid2r& s0, const Grid2r& sTmp, const VelocityField2D<real>& vel, real dt) {
	Vector2<real> pos, posNew;
	real sMax, sMin;

	AABB2D<real> aabb = s.getAABB();

	// Limiting
	for ( index_type i = 1; i < s.nX-1; ++i ) {
		for ( index_type j = 1; j < s.nY-1; ++j ) {
			posNew = TRACER::backtrack(s.pos(i,j), vel, aabb, dt);
			aabb.constrainedBdry(posNew);

			s0.getMaxAndMin(sMax, sMin, posNew);

			if ( sMin < sTmp(i,j) && sTmp(i,j) < sMax ) s(i,j) = sTmp(i,j);	
		}
	}	
}

template<typename TRACER>
void CFDAdvection_SemiLagrangian<TRACER>::revert_sl(Grid3r& s, const Grid3r& s0, const Grid3r& sTmp, const CFDsolver3D& cfd, real dt, MacGrid3D<real>::FIELD_TYPE type) {
	Vector3<real> posNew;
	real sMax, sMin;

	MacGrid3D_VelocityGetter<real> vel(cfd.u, cfd.v, cfd.w, type);

	// Limiting
	for ( index_type i = 1; i < s.nX-1; ++i ) {
		for ( index_type j = 1; j < s.nY-1; ++j ) {
			for ( index_type k = 1; k < s.nZ-1; ++k ) {

				posNew = TRACER::backtrack(s.pos(i,j,k), vel, cfd, dt);
				cfd.constrainedBdry(posNew);

				s0.getMaxAndMin(sMax, sMin, posNew);

				if ( sMin < sTmp(i,j,k) && sTmp(i,j,k) < sMax ) s(i,j,k) = sTmp(i,j,k);			

			}			
		}
	}
}

template<typename TRACER>
void CFDAdvection_SemiLagrangian<TRACER>::revert_sl(Grid3r& s, const Grid3r& s0, const Grid3r& sTmp, const VelocityField3D<real>& vel, real dt) {
	Vector3<real> pos, posNew;
	real sMax, sMin;

	AABB<real> aabb = s.getAABB();

	// Limiting
	for ( index_type i = 1; i < s.nX-1; ++i ) {
		for ( index_type j = 1; j < s.nY-1; ++j ) {
			for ( index_type k = 1; k < s.nZ-1; ++k ) {

				posNew = TRACER::backtrack(s.pos(i,j,k), vel, aabb, dt);
				aabb.constrainedBdry(posNew);

				s0.getMaxAndMin(sMax, sMin, posNew);

				if ( sMin < sTmp(i,j,k) && sTmp(i,j,k) < sMax ) s(i,j,k) = sTmp(i,j,k);	
			}			
		}
	}	
}

};			// fluidspace
