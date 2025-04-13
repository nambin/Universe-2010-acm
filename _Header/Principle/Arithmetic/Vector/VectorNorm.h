#ifndef _ARITHMETIC_VECTOR_NORM_H_
#define _ARITHMETIC_VECTOR_NORM_H_

#include "../../../Basis/Basis/Environments.h"
#include "../../../Basis/ParallelProgramming/ParallelComputable.h"
#include "../Vector/VectorOperationFunctors.h"
#include "../GridLibrary/ArrayFunctor.h"

namespace arithmetic {

template<typename VEC>
class VectorNorm : public basis::ParallelComputable {
public : 
	typedef typename VEC::value_type value_type;
	typedef typename VEC::index_type index_type;

	VectorNorm() {}	
	virtual const value_type operator()(const VEC& x) const = 0;			
};

template<typename VEC>
class VectorNorm_L1 : public VectorNorm<VEC> {
public : 
	VectorNorm_L1() {}
	inline const value_type operator()(const VEC& x) const {
		VectorFunctor_L1Norm<VEC> func(x); 
		compute1D_reduce(func, 0, x.size());		
		return func.getValue();
	}
};

template<typename VEC>
class VectorNorm_L2 : public VectorNorm<VEC> {
public : 
	VectorNorm_L2() {}
	inline const value_type operator()(const VEC& x) const {
		VectorFunctor_L2Norm<VEC> func(x); 
		compute1D_reduce(func, 0, x.size());
		return func.getValue();
	}
};

template<typename VEC>
class VectorNorm_Inf : public VectorNorm<VEC> {
public : 
	VectorNorm_Inf() {}
	inline const value_type operator()(const VEC& x) const {
		VectorFunctor_InfNorm<VEC> func(x); 
		compute1D_reduce(func, 0, x.size());
		return func.getValue();		
	}
};

};		// arithmetic

#endif