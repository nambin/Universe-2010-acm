#ifndef _FLUIDSPACE_CFDSOLVER_3D_H_
#define _FLUIDSPACE_CFDSOLVER_3D_H_

#include "../../../Principle/Geometry/MacGrid/MacGrid3D.h"
#include "../Environment/Environments.h"
#include "../FluidLibrary/VelocityFields.h"
#include "CFDsolver.h"

namespace fluidspace {

using namespace std;

class NATUREDLL CFDsolver3D : public CFDsolver, public geometryspace::MacGrid3D<real>	{
public : 
	typedef geometryspace::MacGrid3D<real>::index_type index_type;
	typedef Grid3<real,real>							grid_type;

public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	CFDsolver3D() {}
	const unsigned int getDimension() const { return 3; }
	const real getCFL(real dt, real h) const;

	const real get_h() const { return h; }

	void update_InternalBoundary(const FluidScene_Objects&);
	void interp_KinematicBoundary(Grid3b&, const MacGrid3D<real>::FIELD_TYPE type) const;

	void setVelocity_InternalBoundary(const FluidScene_Objects*);

	//////////////////////////////////////////////////////////////////////////
	// Initialization
	void set(index_type nX, index_type nY, index_type nZ, real h);
	void setVelocityField(const VelocityField3D<real> &field);

	//////////////////////////////////////////////////////////////////////////
	// IO
	template<typename STREAM> void printMaxVel(STREAM &os = std::cout) const;
	template<typename STREAM> void printMaxDiv(STREAM &os = std::cout) const;

	void printDivergenceInfo(index_type, index_type, index_type) const;

	template<typename GRID3, typename STREAM> void printMaxDiv(const GRID3& grid, typename GRID3::value_type _flag, STREAM &os = cout) const;

public : 
	Grid3b	isKinematicBoundary;
	Grid3r	u,v,w,u0,v0,w0;
	Grid3r  pressure;

public :
	inline real getU(const Vector3<real> &_pos) const { return u.get(_pos);	}
	inline real getV(const Vector3<real> &_pos) const { return v.get(_pos);	}
	inline real getW(const Vector3<real> &_pos) const { return w.get(_pos);	}
	inline Vector3<real> getVel(const Vector3<real> &_pos) const { return Vector3<real>(u.get(_pos), v.get(_pos), w.get(_pos)); }
};

template<typename STREAM>
void CFDsolver3D::printMaxVel(STREAM &os) const {
	os << u.infNorm() << " " << v.infNorm() << " " << w.infNorm() << endl;
}

// Need to Parallelize
template<typename STREAM>
void CFDsolver3D::printMaxDiv(STREAM &os) const {
	struct NegAbsDivergence {
		typedef real ret_type;
		typedef CFDsolver3D::index_type index_type;

		const Grid3r &u, &v, &w;
		const CFDsolver3D &solver;

		NegAbsDivergence(const Grid3r &_u, const Grid3r &_v, const Grid3r &_w, const CFDsolver3D &_s) : u(_u), v(_v), w(_w), solver(_s) {}
		inline real operator()(index_type i, index_type j, index_type k) const { 
			return - ::abs(u(solver.iS2U(i,j,k,1)) - u(solver.iS2U(i,j,k,0)) 
						 + v(solver.iS2V(i,j,k,1)) - v(solver.iS2V(i,j,k,0)) 
						 + w(solver.iS2W(i,j,k,1)) - w(solver.iS2W(i,j,k,0))) * solver.hInv;
		}
	};

	NegAbsDivergence _func(u, v, w, *this);
	TBBFunctor_MinIndex3d<NegAbsDivergence> minIdx(_func);

	minIdx(tbb::blocked_range3d<index_type>(1, nX+1, 1, nY+1, 1, nZ+1));

	char str[512];
	sprintf_s(str, "Maximum Divergence : %30.20g %d %d %d",-NegAbsDivergence(u,v,w,*this)(minIdx.index0, minIdx.index1, minIdx.index2)
														  , minIdx.index0, minIdx.index1, minIdx.index2);
	os << str << std::endl;
}

template<typename GRID3, typename STREAM> 
void CFDsolver3D::printMaxDiv(const GRID3& grid, typename GRID3::value_type _flag, STREAM &os) const {

	real max_div = 0, div = 0;
	index_type ii = 0, jj = 0, kk = 0;

	for(index_type i=1;i<isKinematicBoundary.nX-1;++i) {
		for(index_type j=1;j<isKinematicBoundary.nY-1;++j) {
			for(index_type k=1;k<isKinematicBoundary.nZ-1;++k) {

//				if(grid(i,j,k) == false && isKinematicBoundary(i,j,k) == true)
//					throwError("Error in print div");

				if(grid(i,j,k) != _flag)
					continue;

				div = ::abs(  u(iS2U(i,j,k,1)) - u(iS2U(i,j,k,0)) 
							+ v(iS2V(i,j,k,1)) - v(iS2V(i,j,k,0)) 
							+ w(iS2W(i,j,k,1)) - w(iS2W(i,j,k,0))) * hInv;

				if(max_div < div) {
					max_div = div;
					ii = i; jj = j; kk = k;
				}

			}			
		}
	}

	char str[512];
	sprintf_s(str, "Maximum Divergence for NB       : %30.20g %d %d %d, isKinematicBdry : %d",max_div, ii, jj, kk, isKinematicBoundary(ii,jj,kk));
	os << str << std::endl;

	max_div = 0, div = 0;
	ii = 0, jj = 0, kk = 0;

	for(index_type i=1;i<isKinematicBoundary.nX-1;++i) {
		for(index_type j=1;j<isKinematicBoundary.nY-1;++j) {
			for(index_type k=1;k<isKinematicBoundary.nZ-1;++k) {

				if(grid(i,j,k) != _flag || grid(i-1,j,k) != _flag || grid(i+1,j,k) != _flag
										|| grid(i,j-1,k) != _flag || grid(i,j+1,k) != _flag
										|| grid(i,j,k-1) != _flag || grid(i,j,k+1) != _flag)
					continue;

				div = ::abs(  u(iS2U(i,j,k,1)) - u(iS2U(i,j,k,0)) 
							+ v(iS2V(i,j,k,1)) - v(iS2V(i,j,k,0)) 
							+ w(iS2W(i,j,k,1)) - w(iS2W(i,j,k,0))) * hInv;

				if(max_div < div) {
					max_div = div;
					ii = i; jj = j; kk = k;
				}

			}			
		}
	}

	sprintf_s(str, "Maximum Divergence for NB in NB : %30.20g %d %d %d, isKinematicBdry : %d",max_div, ii, jj, kk, isKinematicBoundary(ii,jj,kk));
	os << str << std::endl;
}

};		// namespace fluidspace

#endif