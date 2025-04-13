#pragma once

#include "SRL2D.h"

namespace geometryspace {

template<typename SRL>
class Upwind_Functor_SRL {
public : 
	typedef typename SRL::index_type	index_type;
	typedef typename SRL::value_type	T;
	typedef typename SRL::value_type	value_type;
	typedef typename SRL				srl_type;

	Upwind_Functor_SRL(const SRL& srl) : s(srl) {}	

protected : 
	const SRL&						s;	
	typename SRL::SRL_Coarse_Lerp	c_lerp;
};

template<typename SRL>
class Upwind_SRL_1st : public Upwind_Functor_SRL<SRL> {
public : 
	Upwind_SRL_1st(const SRL& srl) : Upwind_Functor_SRL(srl) {}

	void dx(T& dx_m, T& dx_p, const SRLNode2D<T>&, unsigned int p[2]) const;
	void dy(T& dy_m, T& dy_p, const SRLNode2D<T>&, unsigned int p[2]) const;

	void dx(T& dx_m, T& dx_p, const SRLNode3D<T>&, unsigned int p[3]) const;
	void dy(T& dy_m, T& dy_p, const SRLNode3D<T>&, unsigned int p[3]) const;
	void dz(T& dz_m, T& dz_p, const SRLNode3D<T>&, unsigned int p[3]) const;

private : 
	mutable T val[1], h[1];		// variable for calculation
};

template<typename SRL>
class Upwind_SRL_2nd : public Upwind_Functor_SRL<SRL> {
public : 
	Upwind_SRL_2nd(const SRL& srl) : Upwind_Functor_SRL(srl) {}

	void dx(T& dx_m, T& dx_p, const SRLNode2D<T>&, unsigned int p[2]) const;
	void dy(T& dy_m, T& dy_p, const SRLNode2D<T>&, unsigned int p[2]) const;

	void dx(T& dx_m, T& dx_p, const SRLNode3D<T>&, unsigned int p[3]) const;
	void dy(T& dy_m, T& dy_p, const SRLNode3D<T>&, unsigned int p[3]) const;
	void dz(T& dz_m, T& dz_p, const SRLNode3D<T>&, unsigned int p[3]) const;

private : 	
	mutable T val[5], h[4];		// variable for calculation
};

template<typename SRL>
class Upwind_SRL_WENO3rd : public Upwind_Functor_SRL<SRL> {
public : 
	Upwind_SRL_WENO3rd(const SRL& srl) : Upwind_Functor_SRL(srl), upwind_1st_srl(srl), upwind_2nd_srl(srl) {}

	void dx(T& dx_m, T& dx_p, const SRLNode2D<T>&, unsigned int p[2]) const;
	void dy(T& dy_m, T& dy_p, const SRLNode2D<T>&, unsigned int p[2]) const;

	void dx(T& dx_m, T& dx_p, const SRLNode3D<T>&, unsigned int p[3]) const;
	void dy(T& dy_m, T& dy_p, const SRLNode3D<T>&, unsigned int p[3]) const;
	void dz(T& dz_m, T& dz_p, const SRLNode3D<T>&, unsigned int p[3]) const;

private : 	
	mutable T val[5], h[4];		// variable for calculation

	Upwind_SRL_1st<SRL> upwind_1st_srl;	
	Upwind_SRL_2nd<SRL> upwind_2nd_srl;
};

template<typename SRL>
class Upwind_SRL_WENO : public Upwind_Functor_SRL<SRL> {
public : 
	Upwind_SRL_WENO(const SRL& srl) : Upwind_Functor_SRL(srl), upwind_1st_srl(srl), upwind_2nd_srl(srl) {}

	void dx(T& dx_m, T& dx_p, const SRLNode2D<T>&, unsigned int p[2]) const;
	void dy(T& dy_m, T& dy_p, const SRLNode2D<T>&, unsigned int p[2]) const;

	void dx(T& dx_m, T& dx_p, const SRLNode3D<T>&, unsigned int p[3]) const;
	void dy(T& dy_m, T& dy_p, const SRLNode3D<T>&, unsigned int p[3]) const;
	void dz(T& dz_m, T& dz_p, const SRLNode3D<T>&, unsigned int p[3]) const;

private : 	
	mutable T val[7], h[6];		// variable for calculation

	Upwind_SRL_1st<SRL> upwind_1st_srl;
	Upwind_SRL_2nd<SRL> upwind_2nd_srl;
};

};		// geometryspace

#include "SRL_UpwindDerivative_Functor_Implementation.h"