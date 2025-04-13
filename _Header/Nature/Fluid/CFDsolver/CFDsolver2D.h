#ifndef _FLUIDSPACE_CFDSOLVER_2D_H_
#define _FLUIDSPACE_CFDSOLVER_2D_H_

#include "../../../Basis/TBBProgramming/TBBFunctions.h"
#include "../../../Principle/Geometry/MacGrid/MacGrid2D.h"
#include "../Environment/Environments.h"
#include "../FluidLibrary/VelocityFields.h"
#include "CFDsolver.h"

namespace fluidspace {

using namespace std;

class NATUREDLL CFDsolver2D : public CFDsolver, public geometryspace::MacGrid2D<real> {
public : 
	typedef geometryspace::MacGrid2D<real>::index_type	index_type;
	typedef Grid2<real,real>							grid_type;

public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	CFDsolver2D() {}
	const unsigned int getDimension() const { return 2; }
	const real getCFL(real dt, real h) const;

	const real get_h() const { return h; }

	void update_InternalBoundary(const FluidScene_Objects&);
	void interp_KinematicBoundary(Grid2b&, const MacGrid2D<real>::FIELD_TYPE type) const;

	void setVelocity_InternalBoundary(const FluidScene_Objects*);

	//////////////////////////////////////////////////////////////////////////
	// Initialization
	void set(index_type nX, index_type nY, real h);
	void setVelocityField(const VelocityField2D<real> &field);

	//////////////////////////////////////////////////////////////////////////
	// Utility
	void setSingleVortexFlow();
	void setRotationFlow(size_t nPerOneRotation);

	//////////////////////////////////////////////////////////////////////////
	// IO
	template<typename STREAM> void printMaxVel(STREAM &os = cout) const;
	template<typename STREAM> void printMaxDiv(STREAM &os = cout) const;

	template<typename GRID2, typename STREAM> void printMaxDiv(const GRID2& grid, typename GRID2::value_type _flag, STREAM &os = cout) const;

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	void glDrawVelocity( real scale = 0.5f, const RGBA &color = RGBA::red(), const RGBA &colorP = RGBA::green(), real pointsize = 1.0f) const;
	void glDrawPressure( const RGBA &color = RGBA::red()) const;
	void glDrawKinematicBoundary( const RGBA &color = RGBA::brown()) const;

	const bool load(const string& filename);
	const bool save(const string& filename) const;

	const std::string name_class() const { return XmlElementName(); }	
	static const std::string XmlElementName() { return std::string("CFDsolver2D"); }

public : 
	Grid2b	isKinematicBoundary;
	Grid2r	u,v,u0,v0;
	Grid2r	pressure;

public :
	inline const real getU(const Vector2<real> &_pos) const { return u.get(_pos); }
	inline const real getV(const Vector2<real> &_pos) const { return v.get(_pos); }
	inline const Vector2<real> getVel(const Vector2<real> &_pos) const { return Vector2<real>(u.get(_pos), v.get(_pos)); }
};

template<typename STREAM> 
void CFDsolver2D::printMaxVel(STREAM &os) const {
	os << "Maximum Velocity : " << u.infNorm() << " " << v.infNorm() << std::endl;
}

template<typename STREAM> 
void CFDsolver2D::printMaxDiv(STREAM &os) const {
	struct NegAbsDivergence {
		typedef real ret_type;
		typedef CFDsolver2D::index_type index_type;

		const Grid2r &u, &v;
		const CFDsolver2D &solver;

		NegAbsDivergence(const Grid2r &_u, const Grid2r &_v, const CFDsolver2D &_s) : u(_u), v(_v), solver(_s) {}
		const real operator()(index_type i, index_type j) const { 
			return - ::abs(u(solver.iS2U(i,j,1)) - u(solver.iS2U(i,j,0)) 
						 + v(solver.iS2V(i,j,1)) - v(solver.iS2V(i,j,0))) * solver.hInv; 
		}
	};

	NegAbsDivergence _func(u, v, *this);
	TBBFunctor_MinIndex2d<NegAbsDivergence> minIdx(_func);

	minIdx(tbb::blocked_range2d<index_type>(1, nX+1, 1, 1, nY+1, 1));

	char str[512];
	sprintf_s(str, "Maximum Divergence : %30.20g %d %d",-NegAbsDivergence(u,v,*this)(minIdx.index0, minIdx.index1), minIdx.index0, minIdx.index1);
	os << str << std::endl;

	//{
	//	int ii = 82, jj = 130;

	//	printf("u+ : %30.20g\n",u[iS2U(ii,jj,1)]);
	//	printf("u- : %30.20g\n",u[iS2U(ii,jj,0)]);
	//	printf("v+ : %30.20g\n",v[iS2V(ii,jj,1)]);
	//	printf("v- : %30.20g\n",v[iS2V(ii,jj,0)]);
	//	printf("div : %30.20g\n", (u[iS2U(ii,jj,1)] - u[iS2U(ii,jj,0)] + v[iS2V(ii,jj,1)] - v[iS2V(ii,jj,0)]) * hInv);	
	//}	
}

template<typename GRID2, typename STREAM> 
void CFDsolver2D::printMaxDiv(const GRID2& grid, typename GRID2::value_type _flag, STREAM &os) const {

	real max_div = 0, div = 0;
	index_type ii = 0, jj = 0;

	for(index_type i=1;i<isKinematicBoundary.nX-1;++i) {
		for(index_type j=1;j<isKinematicBoundary.nY-1;++j) {

			if(grid(i,j) != _flag)
				continue;

			div = ::abs(  u(iS2U(i,j,1)) - u(iS2U(i,j,0)) 
						+ v(iS2V(i,j,1)) - v(iS2V(i,j,0))) * hInv; 

			if(max_div < div) {
				max_div = div;
				ii = i; jj = j;
			}
		}
	}

	char str[512];
	sprintf_s(str, "Maximum Divergence for NB       : %30.20g %d %d", max_div, ii, jj);
	os << str << std::endl;

	max_div = 0, div = 0;
	ii = jj = 0;

	for(index_type i=1;i<isKinematicBoundary.nX-1;++i) {
		for(index_type j=1;j<isKinematicBoundary.nY-1;++j) {

			if(grid(i,j) != _flag || grid(i+1,j) != _flag || grid(i,j+1) != _flag
								  || grid(i-1,j) != _flag || grid(i,j-1) != _flag )
				continue;

			div = ::abs(  u(iS2U(i,j,1)) - u(iS2U(i,j,0)) 
						+ v(iS2V(i,j,1)) - v(iS2V(i,j,0))) * hInv; 

			if(max_div < div) {
				max_div = div;
				ii = i; jj = j;
			}
		}
	}
	
	sprintf_s(str, "Maximum Divergence for NB in NB : %30.20g %d %d", max_div, ii, jj);
	os << str << std::endl;
}

};	// namespace fluidspace

#endif
