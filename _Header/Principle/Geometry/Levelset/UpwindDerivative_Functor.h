#pragma once

#include "../../Arithmetic/GridLibrary/Grid2.h"
#include "../../Arithmetic/GridLibrary/Grid3.h"

namespace geometryspace {

template<typename GRID>
class Upwind_Functor {
public : 
	typedef typename GRID::index_type	index_type;
	typedef typename GRID::value_type	T;
	typedef typename GRID::value_type	value_type;
	typedef typename GRID				grid_type;

	Upwind_Functor(const GRID& grid) : g(grid) {}

	void gather_eno2nd_stencil_x(T v[5], index_type i, index_type j) const {
		v[1] = g(i-1,j); v[2] = g(i,j); v[3] = g(i+1,j);

		v[0] = g(MAX(i-2, 1), j); 
		v[4] = g(MIN(i+2, g.nX-2), j);	
	}
	void gather_eno2nd_stencil_y(T v[5], index_type i, index_type j) const {
		v[1] = g(i,j-1); v[2] = g(i,j); v[3] = g(i,j+1);

		v[0] = g(i, MAX(j-2, 1)); 
		v[4] = g(i, MIN(j+2, g.nY-2));	
	}
	void gather_eno2nd_stencil_x(T v[5], index_type i, index_type j, index_type k) const {
		v[1] = g(i-1,j,k); v[2] = g(i,j,k); v[3] = g(i+1,j,k);

		v[0] = g(MAX(i-2, 1), j,k); 
		v[4] = g(MIN(i+2, g.nX-2), j,k);	
	}
	void gather_eno2nd_stencil_y(T v[5], index_type i, index_type j, index_type k) const {
		v[1] = g(i,j-1,k); v[2] = g(i,j,k); v[3] = g(i,j+1,k);

		v[0] = g(i, MAX(j-2, 1),k); 
		v[4] = g(i, MIN(j+2, g.nY-2),k);	
	}
	void gather_eno2nd_stencil_z(T v[5], index_type i, index_type j, index_type k) const {
		v[1] = g(i,j,k-1); v[2] = g(i,j,k); v[3] = g(i,j,k+1);

		v[0] = g(i, j, MAX(k-2, 1)); 
		v[4] = g(i, j, MIN(k+2, g.nZ-2));	
	}

	void gather_eno_stencil_x(T v[7], index_type i, index_type j) const {
		v[2] = g(i-1,j); v[3] = g(i,j); v[4] = g(i+1,j);

		v[0] = g(MAX(i-3, 1), j); 		v[1] = g(MAX(i-2, 1), j); 
		v[5] = g(MIN(i+2, g.nX-2), j);	v[6] = g(MIN(i+3, g.nX-2), j); 
	}
	void gather_eno_stencil_y(T v[7], index_type i, index_type j) const {
		v[2] = g(i,j-1); v[3] = g(i,j); v[4] = g(i,j+1);

		v[0] = g(i, MAX(j-3, 1)); 		v[1] = g(i, MAX(j-2, 1)); 
		v[5] = g(i, MIN(j+2, g.nY-2));	v[6] = g(i, MIN(j+3, g.nY-2)); 
	}
	void gather_eno_stencil_x(T v[7], index_type i, index_type j, index_type k) const {
		v[2] = g(i-1,j,k); v[3] = g(i,j,k); v[4] = g(i+1,j,k);

		v[0] = g(MAX(i-3, 1), j,k); 		v[1] = g(MAX(i-2, 1), j,k); 
		v[5] = g(MIN(i+2, g.nX-2), j,k);	v[6] = g(MIN(i+3, g.nX-2), j,k); 
	}
	void gather_eno_stencil_y(T v[7], index_type i, index_type j, index_type k) const {
		v[2] = g(i,j-1,k); v[3] = g(i,j,k); v[4] = g(i,j+1,k);

		v[0] = g(i, MAX(j-3, 1),k); 		v[1] = g(i, MAX(j-2, 1),k); 
		v[5] = g(i, MIN(j+2, g.nY-2),k);	v[6] = g(i, MIN(j+3, g.nY-2),k); 
	}

	void gather_eno_stencil_z(T v[7], index_type i, index_type j, index_type k) const {
		v[2] = g(i,j,k-1); v[3] = g(i,j,k); v[4] = g(i,j,k+1);

		v[0] = g(i, j, MAX(k-3, 1)); 		v[1] = g(i, j, MAX(k-2, 1)); 
		v[5] = g(i, j, MIN(k+2, g.nZ-2));	v[6] = g(i, j, MIN(k+3, g.nZ-2)); 
	}

protected : 
	const GRID& g;	
};

template<typename GRID>
class Upwind_1st : public Upwind_Functor<GRID> {
public : 
	Upwind_1st(const GRID& grid) : Upwind_Functor(grid) {}

	const T dx(index_type i, index_type j, bool pos_vel) const { return (pos_vel) ? (g(i,j) - g(i-1,j)) * g.dxInv : (g(i+1,j) - g(i,j)) * g.dxInv; }
	const T dy(index_type i, index_type j, bool pos_vel) const { return (pos_vel) ? (g(i,j) - g(i,j-1)) * g.dyInv : (g(i,j+1) - g(i,j)) * g.dyInv; }

	void dx(T& dx_m, T& dx_p, index_type i, index_type j) const { dx_m = (g(i,j) - g(i-1,j)) * g.dxInv; dx_p = (g(i+1,j) - g(i,j)) * g.dxInv; }
	void dy(T& dy_m, T& dy_p, index_type i, index_type j) const { dy_m = (g(i,j) - g(i,j-1)) * g.dyInv; dy_p = (g(i,j+1) - g(i,j)) * g.dyInv; }

	const T dx(index_type i, index_type j, index_type k, bool pos_vel) const { return (pos_vel) ? (g(i,j,k) - g(i-1,j,k)) * g.dxInv : (g(i+1,j,k) - g(i,j,k)) * g.dxInv; }
	const T dy(index_type i, index_type j, index_type k, bool pos_vel) const { return (pos_vel) ? (g(i,j,k) - g(i,j-1,k)) * g.dyInv : (g(i,j+1,k) - g(i,j,k)) * g.dyInv; }
	const T dz(index_type i, index_type j, index_type k, bool pos_vel) const { return (pos_vel) ? (g(i,j,k) - g(i,j,k-1)) * g.dzInv : (g(i,j,k+1) - g(i,j,k)) * g.dzInv; }

	void dx(T& dx_m, T& dx_p, index_type i, index_type j, index_type k) const { dx_m = (g(i,j,k) - g(i-1,j,k)) * g.dxInv; dx_p = (g(i+1,j,k) - g(i,j,k)) * g.dxInv; }
	void dy(T& dy_m, T& dy_p, index_type i, index_type j, index_type k) const { dy_m = (g(i,j,k) - g(i,j-1,k)) * g.dyInv; dy_p = (g(i,j+1,k) - g(i,j,k)) * g.dyInv; }
	void dz(T& dz_m, T& dz_p, index_type i, index_type j, index_type k) const { dz_m = (g(i,j,k) - g(i,j,k-1)) * g.dzInv; dz_p = (g(i,j,k+1) - g(i,j,k)) * g.dzInv; }
};

template<typename GRID>
class Upwind_2nd : public Upwind_Functor<GRID> {
public : 
	Upwind_2nd(const GRID& grid) : Upwind_Functor(grid) {}

	const T dx(index_type i, index_type j, bool pos_vel) const { gather_eno2nd_stencil_x(v, i, j); return eno2nd(v, g.dx, pos_vel); }
	const T dy(index_type i, index_type j, bool pos_vel) const { gather_eno2nd_stencil_y(v, i, j); return eno2nd(v, g.dy, pos_vel);	}

	void dx(T& dx_m, T& dx_p, index_type i, index_type j) const { gather_eno2nd_stencil_x(v, i, j);	dx_m = eno2nd(v, g.dx, true); dx_p = eno2nd(v, g.dx, false); }
	void dy(T& dy_m, T& dy_p, index_type i, index_type j) const { gather_eno2nd_stencil_y(v, i, j); dy_m = eno2nd(v, g.dy, true); dy_p = eno2nd(v, g.dy, false); }

	const T dx(index_type i, index_type j, index_type k, bool pos_vel) const { gather_eno2nd_stencil_x(v, i, j, k); return eno2nd(v, g.dx, pos_vel); }
	const T dy(index_type i, index_type j, index_type k, bool pos_vel) const { gather_eno2nd_stencil_y(v, i, j, k); return eno2nd(v, g.dy, pos_vel); }
	const T dz(index_type i, index_type j, index_type k, bool pos_vel) const { gather_eno2nd_stencil_z(v, i, j, k); return eno2nd(v, g.dz, pos_vel); }

	void dx(T& dx_m, T& dx_p, index_type i, index_type j, index_type k) const {	gather_eno2nd_stencil_x(v, i, j, k); dx_m = eno2nd(v, g.dx, true); dx_p = eno2nd(v, g.dx, false); }
	void dy(T& dy_m, T& dy_p, index_type i, index_type j, index_type k) const { gather_eno2nd_stencil_y(v, i, j, k); dy_m = eno2nd(v, g.dy, true); dy_p = eno2nd(v, g.dy, false); }
	void dz(T& dz_m, T& dz_p, index_type i, index_type j, index_type k) const { gather_eno2nd_stencil_z(v, i, j, k); dz_m = eno2nd(v, g.dz, true); dz_p = eno2nd(v, g.dz, false); }

private : 
	mutable T v[7];
};

template<typename GRID>
class Upwind_ENO : public Upwind_Functor<GRID> {
public : 
	Upwind_ENO(const GRID& grid) : Upwind_Functor(grid) {}

	const T dx(index_type i, index_type j, bool pos_vel) const { gather_eno_stencil_x(v, i, j); return eno(v, g.dx, pos_vel); }
	const T dy(index_type i, index_type j, bool pos_vel) const { gather_eno_stencil_y(v, i, j); return eno(v, g.dy, pos_vel); }

	void dx(T& dx_m, T& dx_p, index_type i, index_type j) const { gather_eno_stencil_x(v, i, j); dx_m = eno(v, g.dx, true); dx_p = eno(v, g.dx, false); }
	void dy(T& dy_m, T& dy_p, index_type i, index_type j) const { gather_eno_stencil_y(v, i, j); dy_m = eno(v, g.dy, true);	dy_p = eno(v, g.dy, false);	}

	const T dx(index_type i, index_type j, index_type k, bool pos_vel) const { gather_eno_stencil_x(v, i, j, k); return eno(v, g.dx, pos_vel); }
	const T dy(index_type i, index_type j, index_type k, bool pos_vel) const { gather_eno_stencil_y(v, i, j, k); return eno(v, g.dy, pos_vel); }
	const T dz(index_type i, index_type j, index_type k, bool pos_vel) const { gather_eno_stencil_z(v, i, j, k); return eno(v, g.dz, pos_vel); }

	void dx(T& dx_m, T& dx_p, index_type i, index_type j, index_type k) const {	gather_eno_stencil_x(v, i, j, k); dx_m = eno(v, g.dx, true); dx_p = eno(v, g.dx, false); }
	void dy(T& dy_m, T& dy_p, index_type i, index_type j, index_type k) const { gather_eno_stencil_y(v, i, j, k); dy_m = eno(v, g.dy, true); dy_p = eno(v, g.dy, false); }
	void dz(T& dz_m, T& dz_p, index_type i, index_type j, index_type k) const { gather_eno_stencil_z(v, i, j, k); dz_m = eno(v, g.dz, true); dz_p = eno(v, g.dz, false); }

private : 
	mutable T v[7];
};

template<typename GRID>
class Upwind_WENO : public Upwind_Functor<GRID> {
public : 
	Upwind_WENO(const GRID& grid) : Upwind_Functor(grid) {}

	const T dx(index_type i, index_type j, bool pos_vel) const { gather_eno_stencil_x(v, i, j); return weno(v, g.dx, pos_vel); }
	const T dy(index_type i, index_type j, bool pos_vel) const { gather_eno_stencil_y(v, i, j); return weno(v, g.dy, pos_vel); }

	void dx(T& dx_m, T& dx_p, index_type i, index_type j) const { gather_eno_stencil_x(v, i, j); dx_m = weno(v, g.dx, true); dx_p = weno(v, g.dx, false); }
	void dy(T& dy_m, T& dy_p, index_type i, index_type j) const { gather_eno_stencil_y(v, i, j); dy_m = weno(v, g.dy, true); dy_p = weno(v, g.dy, false); }

	const T dx(index_type i, index_type j, index_type k, bool pos_vel) const { gather_eno_stencil_x(v, i, j, k); return weno(v, g.dx, pos_vel); }
	const T dy(index_type i, index_type j, index_type k, bool pos_vel) const { gather_eno_stencil_y(v, i, j, k); return weno(v, g.dy, pos_vel); }
	const T dz(index_type i, index_type j, index_type k, bool pos_vel) const { gather_eno_stencil_z(v, i, j, k); return weno(v, g.dz, pos_vel); }

	void dx(T& dx_m, T& dx_p, index_type i, index_type j, index_type k) const {	gather_eno_stencil_x(v, i, j, k); dx_m = weno(v, g.dx, true); dx_p = weno(v, g.dx, false); }
	void dy(T& dy_m, T& dy_p, index_type i, index_type j, index_type k) const { gather_eno_stencil_y(v, i, j, k); dy_m = weno(v, g.dy, true); dy_p = weno(v, g.dy, false); }
	void dz(T& dz_m, T& dz_p, index_type i, index_type j, index_type k) const { gather_eno_stencil_z(v, i, j, k); dz_m = weno(v, g.dz, true); dz_p = weno(v, g.dz, false); }

private : 
	mutable T v[7];
};

};		// geometryspace