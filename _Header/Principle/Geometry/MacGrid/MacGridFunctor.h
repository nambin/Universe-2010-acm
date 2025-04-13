#ifndef _GEOMETRYSPACE_MACGRID_FUNCTOR_H_
#define _GEOMETRYSPACE_MACGRID_FUNCTOR_H_

namespace geometryspace {

template<typename T>
class MacGrid2D_VelocityGetter {
public : 
	typedef typename MacGrid2D<T>::index_type index_type;
	
public : 
	const Grid2<T,T>					& u, & v;
	typename MacGrid2D<T>::FIELD_TYPE	type;

	MacGrid2D<T>						mac;

	MacGrid2D_VelocityGetter(const Grid2<T,T> &_u, const Grid2<T,T> &_v, typename MacGrid2D<T>::FIELD_TYPE _type)
	: u(_u), v(_v), type(_type), mac(u.nX-1, u.nY-2, u.dx) { assert(u.nX+1==v.nX); assert(u.nY==v.nY+1); }

	const Vector2<T> getVel(index_type i, index_type j) const {
		switch(type) {
			case MacGrid2D<T>::TYPE_U :	return Vector2<T>(u(i,j), 0.25f*(v(i,j-1)+v(i,j)+v(i+1,j-1)+v(i+1,j)));
										break;
			case MacGrid2D<T>::TYPE_V :	return Vector2<T>(0.25f*(u(i,j+1)+u(i,j)+u(i-1,j+1)+u(i-1,j)), v(i,j));
										break;
			case MacGrid2D<T>::TYPE_S :	return Vector2<T>((u(mac.iS2U(i,j,0)) + u(mac.iS2U(i,j,1))) * 0.5f
													    , (v(mac.iS2V(i,j,0)) + v(mac.iS2V(i,j,1))) * 0.5f);
										break;
		}
		return Vector2<T>(0);
	}
	const Vector2<T> getVel(const Vector2<T>& pos) const { return Vector2<T>(u(pos), v(pos)); }
};

template<typename T>
class MacGrid2D_advect_SemiLagrangian {
	Grid2<T,T>							& s;
	const Grid2<T,T>					& s0;
	T									dt; 
	
	const MacGrid2D<T>					& macGrid;

	const MacGrid2D_VelocityGetter<T>	vel;

public : 
	typedef typename MacGrid2D<T>::index_type index_type;

public : 
	MacGrid2D_advect_SemiLagrangian(Grid2<T,T> &_s, const Grid2<T,T> &_s0, const Grid2<T,T> &_u_Base, const Grid2<T,T> &_v_Base
								  , T _dt, typename MacGrid2D<T>::FIELD_TYPE _type, const MacGrid2D<T> &_mac)
	: s(_s), s0(_s0), dt(_dt), macGrid(_mac), vel(_u_Base, _v_Base, _type) {}

	void operator()(const basis::blocked_range2d<index_type>& r) const {

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
				posNew = s.pos(i,j) - dt * vel.getVel(i,j);
				macGrid.constrainedBdry(posNew);

				s(i,j) = s0.get(posNew);
			}
		}
	}
};

template<typename T>
class MacGrid3D_VelocityGetter {
public : 
	typedef typename MacGrid3D<T>::index_type index_type;

public : 
	const Grid3<T,T>					& u, & v, & w;
	typename MacGrid3D<T>::FIELD_TYPE	type;

	MacGrid3D<T>						mac;

	MacGrid3D_VelocityGetter(const Grid3<T,T> &_u, const Grid3<T,T> &_v, const Grid3<T,T> &_w, typename MacGrid3D<T>::FIELD_TYPE _type)
		: u(_u), v(_v), w(_w), type(_type), mac(u.nX-1, u.nY-2, u.nZ-2, u.dx) { assert(u.nX+1==v.nX); assert(u.nY==v.nY+1); assert(u.nZ==v.nZ); }

	const Vector3<T> getVel(index_type i, index_type j, index_type k) const {
		T un, vn, wn;
		switch(type) {
			case MacGrid3D<T>::TYPE_U :	un = u(i,j,k);
										vn = T(0.25)*(v(i,j-1,k)+v(i,j,k)+v(i+1,j-1,k)+v(i+1,j,k));
										wn = T(0.25)*(w(i,j,k-1)+w(i,j,k)+w(i+1,j,k-1)+w(i+1,j,k));
										break;
			case MacGrid3D<T>::TYPE_V :	un = T(0.25)*(u(i-1,j,k)+u(i,j,k)+u(i-1,j+1,k)+u(i,j+1,k));
										vn = v(i,j,k);
										wn = T(0.25)*(w(i,j,k-1)+w(i,j,k)+w(i,j+1,k-1)+w(i,j+1,k));
										break;
			case MacGrid3D<T>::TYPE_W :	un = T(0.25)*(u(i-1,j,k)+u(i,j,k)+u(i-1,j,k+1)+u(i,j,k+1));
										vn = T(0.25)*(v(i,j-1,k)+v(i,j,k)+v(i,j-1,k+1)+v(i,j,k+1));
										wn = w(i,j,k);
										break;
			case MacGrid3D<T>::TYPE_S :	un = T(0.5) * (u[mac.iS2U(i,j,k,0)] + u[mac.iS2U(i,j,k,1)]);
										vn = T(0.5) * (v[mac.iS2V(i,j,k,0)] + v[mac.iS2V(i,j,k,1)]);
										wn = T(0.5) * (w[mac.iS2W(i,j,k,0)] + w[mac.iS2W(i,j,k,1)]);
										break;
		}
		return Vector3<T>(un,vn,wn);
	}
	const Vector3<T> getVel(const Vector3<T>& pos) const { return Vector3<T>(u(pos), v(pos), w(pos)); }
};

template<typename T>
class MacGrid3D_advect_SemiLagrangian {
	Grid3<T,T>							& s;
	const Grid3<T,T>					& s0;
	const MacGrid3D<T>					& macGrid;

	T									dt; 	
	const MacGrid3D_VelocityGetter<T>	vel;

public : 
	typedef typename MacGrid3D<T>::index_type index_type;

public : 
	MacGrid3D_advect_SemiLagrangian(Grid3<T,T> &_s, const Grid3<T,T> &_s0, const Grid3<T,T> &_u_Base, const Grid3<T,T> &_v_Base, const Grid3<T,T> &_w_Base
		, T _dt, typename MacGrid3D<T>::FIELD_TYPE _type, const MacGrid3D<T> &_mac)
		: s(_s), s0(_s0), vel(_u_Base, _v_Base, _w_Base, _type), dt(_dt), macGrid(_mac) {}

	void operator()( const basis::blocked_range3d<index_type>& r ) const {	

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
					posNew = s.pos(i,j,k) - dt * vel.getVel(i,j,k);
					macGrid.constrainedBdry(posNew);
					s(i,j,k) = s0.get(posNew);
				}
			}
		}
	}
};

};		// geometryspace

#endif