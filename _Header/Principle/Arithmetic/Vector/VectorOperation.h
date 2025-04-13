#ifndef _ARITHMETIC_VECTOR_OPERATION_H_
#define _ARITHMETIC_VECTOR_OPERATION_H_

#include "../../../Basis/Basis/Environments.h"
#include "../../../Basis/ParallelProgramming/ParallelComputable.h"
#include "VectorOperationFunctors.h"
#include "VectorNorm.h"

namespace arithmetic {

#ifdef _DEBUG
#define __ASSERT(c) assert(c)
#else
#define __ASSERT(c) {}
#endif

template<typename VEC>
class VectorOperation : public basis::ParallelComputable {
public : 
	typedef typename VEC::index_type index_type;

public :
	VectorOperation() : basis::ParallelComputable() {}
	VectorOperation(const basis::ParallelComputable& c) : basis::ParallelComputable(c) {}
	
	inline void zero(VEC& x) const {	
		set(x, typename VEC::value_type(0.0));
	}

	inline void set(VEC& x, typename const VEC::value_type val) const {	
		compute1D_for(VectorFunctor_SetValue<VEC, typename VEC::value_type>(x,val), 0, x.size());
	}	

	template<typename VECA>
	inline void set(VEC& x, const VECA& a) const {
		compute1D_for(VectorFunctor_Set<VEC,VECA>(x,a), 0, x.size());
	}

	template<typename VECA>
	inline void neg(VEC& x, const VECA& a) const {
		__ASSERT( x.size() == a.size() );	
		compute1D_for(VectorFunctor_NegativeSet<VEC,VECA>(x,a), 0, x.size());
	}	

	template<typename VECA, typename VECB>
	inline void add(VEC& x, const VECA& a, const VECB& b) const {
		__ASSERT( x.size() == a.size() && a.size() == b.size() );	
		compute1D_for(VectorFunctor_Add<VEC,VECA,VECB>(x,a,b), 0, x.size());
	}	

	template<typename VECA, typename VECB>
	inline void dif(VEC& x, const VECA& a, const VECB& b) const {
		__ASSERT( x.size() == a.size() && a.size() == b.size() );
		compute1D_for(VectorFunctor_Dif<VEC,VECA,VECB>(x,a,b), 0, x.size());
	}

	template<typename VECA>
	inline void icr(VEC& x, const VECA& a) const {
		__ASSERT( x.size() == a.size() );	
		compute1D_for(VectorFunctor_Icr<VEC,VECA>(x,a), 0, x.size());		
	}

	template<typename VECA>
	inline void dcr(VEC& x, const VECA& a) const {
		__ASSERT( x.size() == a.size() );	
		compute1D_for(VectorFunctor_Dcr<VEC,VECA>(x,a), 0, x.size());		
	}

	template<typename VECA, typename R>
	inline void mul(VEC&x, const VECA& a, R c) const {
		__ASSERT( x.size() == a.size() );	
		compute1D_for(VectorFunctor_Mul<VEC,VECA,R>(x,a,c), 0, x.size());	
	}

	template<typename R>
	inline void scale(VEC&x, R c) const {
		compute1D_for(VectorFunctor_Scale<VEC,R>(x,c), 0, x.size());	
	}

	template<typename VECY>
	inline const typename VEC::value_type dot(const VEC& x, const VECY& y) const {
		__ASSERT( x.size() == y.size() );	
		VectorFunctor_Dot<VEC,VECY> func(x,y); 
		compute1D_reduce(func, 0, x.size());		
		return func.getValue();
	}	

	inline const typename VEC::value_type dot_directly(const VEC& x, const VEC& y) const {
		__ASSERT( x.size() == y.size() );	
		VectorFunctor_Dot_Directly<VEC> func(x,y); 
		compute1D_reduce(func, 0, x.size());		
		return func.getValue();
	}	

	template<typename VECA, typename R>
	inline void icr_scaled(VEC&x, const VECA& a, R c) const {
		__ASSERT( x.size() == a.size() );	
		compute1D_for(VectorFunctor_IcrScaled<VEC,VECA,R>(x,a,c), 0, x.size());
	}

	inline void icr_scaled_directly(VEC&x, const VEC& a, typename VEC::value_type c) const {
		__ASSERT( x.size() == a.size() );	
		compute1D_for(VectorFunctor_IcrScaled_Directly<VEC>(x,a,c), 0, x.size());
	}

	template<typename VECA, typename R>
	inline void scaled_icr(VEC&x, R c, const VECA& a) const {
		__ASSERT( x.size() == a.size() );	
		compute1D_for(VectorFunctor_ScaledIcr<VEC,VECA,R>(x,a,c), 0, x.size());
	}

	inline void scaled_icr_directly(VEC&x, typename VEC::value_type c, const VEC& a) const {
		__ASSERT( x.size() == a.size() );	
		compute1D_for(VectorFunctor_ScaledIcr_Directly<VEC>(x,a,c), 0, x.size());
	}

	inline const typename VEC::value_type sum(const VEC& x) const {
		VectorFunctor_Sum<VEC> func(x); 
		compute1D_reduce(func, 0, x.size());
		return func.getValue();		
	}

	template<typename VECY>
	inline const typename VEC::value_type length(const VEC& a, const VECY& b) const {
		return sqrt(sqrlength(a,b));
	}

	template<typename VECY>
	inline const typename VEC::value_type sqrlength(const VEC& x, const VECY& y) const {
		__ASSERT( x.size() == y.size() );	
		VectorFunctor_SqrLength<VEC,VECY> func(x,y); 
		compute1D_reduce(func, 0, x.size());
		return func.getValue();			
	}

	template<typename VECA, typename VECB>
	inline void maxVec(VEC& x, const VECA& a, const VECB& b) const {
		__ASSERT( x.size() == a.size() && x.size() == b.size() );		
		compute1D_for(VectorFunctor_MaxVec<VEC,VECA,VECB>(x,a,b), 0, x.size());
	}

	template<typename VECA, typename VECB>
	inline void minVec(VEC& x, const VECA& a, const VECB& b) const {
		__ASSERT( x.size() == a.size() && x.size() == b.size() );		
		compute1D_for(VectorFunctor_MinVec<VEC,VECA,VECB>(x,a,b), 0, x.size());
	}

	inline const typename VEC::value_type l1Norm(const VEC& v) const { return VectorNorm_L1<VEC>(*this).norm(v);  }
	inline const typename VEC::value_type l2Norm(const VEC& v) const { return VectorNorm_L2<VEC>(*this).norm(v);  }
	inline const typename VEC::value_type infNorm(const VEC& v) const { return VectorNorm_Inf<VEC>(*this).norm(v); }
};

// -V
template <typename VECTOR>
struct Vector_neg_V {
	VECTOR const & V ;
	Vector_neg_V(VECTOR const & _V) : V(_V) {}
} ;

// V + U
template <typename VECTOR>
struct Vector_V_add_V {
	VECTOR const & V1 ;
	VECTOR const & V2 ;
	Vector_V_add_V(VECTOR const & _V1, VECTOR const & _V2) : V1(_V1), V2(_V2) {}
} ;

// V - U
template <typename VECTOR>
struct Vector_V_dif_V {
	VECTOR const & V1 ;
	VECTOR const & V2 ;
	Vector_V_dif_V(VECTOR const & _V1, VECTOR const & _V2) : V1(_V1), V2(_V2) {}
} ;

// V * s
template <typename VECTOR, typename SCALAR>
struct Vector_V_mul_S {
	VECTOR const & V ;
	SCALAR const & s ;
	Vector_V_mul_S(VECTOR const & _V, SCALAR const & _s) : V(_V), s(_s) {}
} ;

#undef __ASSERT

};		// arithmetic

#endif

