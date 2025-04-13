#pragma once

namespace fluidspace {

template<typename T, typename VELFIELD, typename TRACER> 
class VelocityField2D_advect_SemiLagrangian {
	Grid2<T,T>							& s;
	const Grid2<T,T>					& s0;
	const VELFIELD						& vel;
	T									dt; 

public : 
	typedef typename Grid2<T,T>::index_type index_type;

	VelocityField2D_advect_SemiLagrangian(Grid2<T,T> &_s, const Grid2<T,T> &_s0, const VELFIELD& _vel, T _dt)
	: s(_s), s0(_s0), vel(_vel), dt(_dt) {}

	void operator()(const basis::blocked_range2d<index_type>& r) const {		
		const AABB2D<T> aabb = s.getAABB();		

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
			Vector2<T> posNew;

			for( index_type j=r.cols().begin(); j!=r.cols().end( ); ++j ) {
				posNew = TRACER::backtrack(s.pos(i,j), vel, aabb, dt);
				aabb.constrainedBdry(posNew);

				s(i,j) = s0.get(posNew);
			}
		}
	}
};

template<typename T, typename VELFIELD, typename TRACER> 
class VelocityField3D_advect_SemiLagrangian {
	Grid3<T,T>							& s;
	const Grid3<T,T>					& s0;
	const VELFIELD						& vel;
	T									dt; 

public : 
	typedef typename Grid3<T,T>::index_type index_type;

	VelocityField3D_advect_SemiLagrangian(Grid3<T,T> &_s, const Grid3<T,T> &_s0, const VELFIELD& _vel, T _dt)
	: s(_s), s0(_s0), vel(_vel), dt(_dt) {}

	void operator()(const basis::blocked_range3d<index_type>& r) const {		
		const AABB<T> aabb = s.getAABB();

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

			for( index_type j=r.rows().begin(); j!=r.rows().end(); ++j ) {
				for( index_type k=r.cols().begin(); k!=r.cols().end(); ++k ) {
					posNew = TRACER::backtrack(s.pos(i,j,k), vel, aabb, dt);
					aabb.constrainedBdry(posNew);

					s(i,j,k) = s0.get(posNew);
				}
			}
		}
	}
};

};		// fluidspace