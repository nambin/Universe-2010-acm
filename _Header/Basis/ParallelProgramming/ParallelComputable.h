#ifndef _BASIS_PARALLEL_COMPUTABLE_H_
#define _BASIS_PARALLEL_COMPUTABLE_H_

#include "../Basis.h"
#include "../Basis/Environments.h"
#include "../Basis/blocked_range.h"
#include "../Basis/blocked_range2d.h"
#include "../Basis/blocked_range3d.h"

namespace basis {

class BASISDLL ParallelComputable {
public : 
	ParallelComputable() {}
	ParallelComputable(const ParallelComputable& c) { setComputingOption(c); }

	virtual void setComputingOption(const ParallelComputable& rhs) { (*this) = rhs; }

	template<typename FUNCTOR> inline void compute1D_for(FUNCTOR& func, typename FUNCTOR::index_type begin, typename FUNCTOR::index_type end) const;
	template<typename FUNCTOR> inline void compute2D_for(FUNCTOR& func, typename FUNCTOR::index_type bx, typename FUNCTOR::index_type ex, typename FUNCTOR::index_type by, typename FUNCTOR::index_type ey) const;
	template<typename FUNCTOR> inline void compute3D_for(FUNCTOR& func, typename FUNCTOR::index_type bx, typename FUNCTOR::index_type ex, typename FUNCTOR::index_type by, typename FUNCTOR::index_type ey, typename FUNCTOR::index_type bz, typename FUNCTOR::index_type ez) const;

	template<typename FUNCTOR> inline void compute1D_reduce(FUNCTOR& func, typename FUNCTOR::index_type begin, typename FUNCTOR::index_type end) const;
	template<typename FUNCTOR> inline void compute2D_reduce(FUNCTOR& func, typename FUNCTOR::index_type bx, typename FUNCTOR::index_type ex, typename FUNCTOR::index_type by, typename FUNCTOR::index_type ey) const;
	template<typename FUNCTOR> inline void compute3D_reduce(FUNCTOR& func, typename FUNCTOR::index_type bx, typename FUNCTOR::index_type ex, typename FUNCTOR::index_type by, typename FUNCTOR::index_type ey, typename FUNCTOR::index_type bz, typename FUNCTOR::index_type ez) const;
};

template<typename FUNCTOR>
inline void ParallelComputable::compute1D_for(FUNCTOR& func, typename FUNCTOR::index_type begin, typename FUNCTOR::index_type end) const {
	func(blocked_range<typename FUNCTOR::index_type>(begin,end));
}

template<typename FUNCTOR>
inline void ParallelComputable::compute2D_for(FUNCTOR& func, typename FUNCTOR::index_type bx, typename FUNCTOR::index_type ex, typename FUNCTOR::index_type by, typename FUNCTOR::index_type ey) const {
	func(blocked_range2d<typename FUNCTOR::index_type>(bx,ex,by,ey));
}

template<typename FUNCTOR>
inline void ParallelComputable::compute3D_for(FUNCTOR& func, typename FUNCTOR::index_type bx, typename FUNCTOR::index_type ex, typename FUNCTOR::index_type by, typename FUNCTOR::index_type ey, typename FUNCTOR::index_type bz, typename FUNCTOR::index_type ez) const {
	func(blocked_range3d<typename FUNCTOR::index_type>(bx,ex,by,ey,bz,ez));
}

template<typename FUNCTOR>
inline void ParallelComputable::compute1D_reduce(FUNCTOR& func, typename FUNCTOR::index_type begin, typename FUNCTOR::index_type end) const {
	func(blocked_range<typename FUNCTOR::index_type>(begin,end));
}

template<typename FUNCTOR>
inline void ParallelComputable::compute2D_reduce(FUNCTOR& func, typename FUNCTOR::index_type bx, typename FUNCTOR::index_type ex, typename FUNCTOR::index_type by, typename FUNCTOR::index_type ey) const {
	func(blocked_range2d<typename FUNCTOR::index_type>(bx,ex,by,ey));
}

template<typename FUNCTOR>
inline void ParallelComputable::compute3D_reduce(FUNCTOR& func, typename FUNCTOR::index_type bx, typename FUNCTOR::index_type ex, typename FUNCTOR::index_type by, typename FUNCTOR::index_type ey, typename FUNCTOR::index_type bz, typename FUNCTOR::index_type ez) const {
	func(blocked_range3d<typename FUNCTOR::index_type>(bx,ex,by,ey,bz,ez));
}

};		// basis

#endif