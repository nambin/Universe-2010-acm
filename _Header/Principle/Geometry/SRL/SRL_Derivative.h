#pragma once

#include "SRL2D.h"

namespace geometryspace {

template<typename SRL>
class Derivative_Functor_SRL {
public : 
	typedef typename SRL::index_type	index_type;
	typedef typename SRL::value_type	T;
	typedef typename SRL::value_type	value_type;
	typedef typename SRL				srl_type;

	Derivative_Functor_SRL(const SRL& srl) : s(srl) {}	

protected : 
	const SRL&						s;	
	typename SRL::SRL_Coarse_Lerp	c_lerp;
};

template<typename SRL>
class Derivative_SRL_2nd : public Derivative_Functor_SRL<SRL> {
public : 
	Derivative_SRL_2nd(const SRL& srl) : Derivative_Functor_SRL(srl) {}

	const T dx(const SRLNode3D<T>&, unsigned int p[3]) const;
	const T dy(const SRLNode3D<T>&, unsigned int p[3]) const;
	const T dz(const SRLNode3D<T>&, unsigned int p[3]) const;

private : 	
	mutable T val[3], h[2];		// variable for calculation
};

};		// geometryspace

#include "SRL_Derivative_Functor.h"