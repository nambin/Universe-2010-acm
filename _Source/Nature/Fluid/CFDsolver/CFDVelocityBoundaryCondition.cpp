#include <stdafx.h>
#include <Nature/Fluid/CFDsolver/CFDVelocityBoundaryCondition.h>
#include <Basis/TBBProgramming/TBBFunctions.h>

using namespace std;
using namespace arithmetic;
using namespace geometryspace;
using namespace fluidspace;
using namespace basis;

void CFDVelocityBoundaryCondition::setComputingOption(const ParallelComputable& rhs) {
	ParallelComputable::setComputingOption(rhs);
	if(uBdry)	uBdry->setComputingOption(rhs);
	if(vBdry)	vBdry->setComputingOption(rhs);
	if(wBdry)	wBdry->setComputingOption(rhs);
}

void CFDVelocityBoundaryCondition::enforce(CFDsolver2D& cfd) {
	if(uBdry)	uBdry->enforce(cfd.u);
	if(vBdry)	vBdry->enforce(cfd.v);
}

void CFDVelocityBoundaryCondition::enforce(CFDsolver3D& cfd) {
	if(uBdry)	uBdry->enforce(cfd.u);
	if(vBdry)	vBdry->enforce(cfd.v);
	if(wBdry)	wBdry->enforce(cfd.w);
}

void GridBoundaryCondition_U_NoSlip::enforce(Grid2<real,real>& u) const {
	compute1D_for(Grid2_XDirichlet<real,real>(u, 0), 0, u.nY);
	compute1D_for(Grid2_YNeumann_Scaling<real,real>(u,-1), 0, u.nX);
}

void GridBoundaryCondition_U_NoSlip::enforce(Grid3<real,real>& u) const {
	Grid3_XDirichlet<real, real>(u, 0)(basis::blocked_range2d<grid_index_type>(0, u.nY, 0, u.nZ));
	Grid3_YNeumann_Scaling<real, real>(u,-1)(basis::blocked_range2d<grid_index_type>(0, u.nX, 0, u.nZ));
	Grid3_ZNeumann_Scaling<real, real>(u,-1)(basis::blocked_range2d<grid_index_type>(0, u.nX, 0, u.nY));
}

void GridBoundaryCondition_V_NoSlip::enforce(Grid2<real,real>& v) const {
	compute1D_for(Grid2_YDirichlet<real,real>(v, 0), 0, v.nX);
	compute1D_for(Grid2_XNeumann_Scaling<real,real>(v,-1), 0, v.nY);
}

void GridBoundaryCondition_V_NoSlip::enforce(Grid3<real,real>& v) const {
	Grid3_YDirichlet<real, real>(v, 0)(blocked_range2d<grid_index_type>(0, v.nX, 0, v.nZ));
	Grid3_XNeumann_Scaling<real, real>(v,-1)(blocked_range2d<grid_index_type>(0, v.nY, 0, v.nZ));
	Grid3_ZNeumann_Scaling<real, real>(v,-1)(blocked_range2d<grid_index_type>(0, v.nX, 0, v.nY));
}

void GridBoundaryCondition_W_NoSlip::enforce(Grid3<real,real>& w) const {
	Grid3_ZDirichlet<real, real>(w, 0)(blocked_range2d<grid_index_type>(0, w.nX, 0, w.nY));
	Grid3_XNeumann_Scaling<real, real>(w,-1)(blocked_range2d<grid_index_type>(0, w.nY, 0, w.nZ));
	Grid3_YNeumann_Scaling<real, real>(w,-1)(blocked_range2d<grid_index_type>(0, w.nX, 0, w.nZ));
}

void GridBoundaryCondition_U_FreeSlip::enforce(Grid2<real,real>& u) const {
	compute1D_for(Grid2_XDirichlet<real,real>(u, 0), 0, u.nY);
	compute1D_for(Grid2_YNeumann<real,real>(u), 0, u.nX);	
}

void GridBoundaryCondition_U_FreeSlip::enforce(Grid3<real,real>& u) const {
	Grid3_XDirichlet<real, real>(u, 0)(blocked_range2d<grid_index_type>(0, u.nY, 0, u.nZ));
	Grid3_YNeumann<real, real> uy(u);	uy(blocked_range2d<grid_index_type>(0, u.nX, 0, u.nZ));
	Grid3_ZNeumann<real, real> uz(u);	uz(blocked_range2d<grid_index_type>(0, u.nX, 0, u.nY));
}

void GridBoundaryCondition_V_FreeSlip::enforce(Grid2<real,real>& v) const {
	compute1D_for(Grid2_YDirichlet<real,real>(v, 0), 0, v.nX);
	compute1D_for(Grid2_XNeumann<real,real>(v), 0, v.nY);
}

void GridBoundaryCondition_V_FreeSlip::enforce(Grid3<real,real>& v) const {
	Grid3_YDirichlet<real, real>(v, 0)(blocked_range2d<grid_index_type>(0, v.nX, 0, v.nZ));
	Grid3_XNeumann<real, real> vx(v);	vx(blocked_range2d<grid_index_type>(0, v.nY, 0, v.nZ));
	Grid3_ZNeumann<real, real> vz(v);	vz(blocked_range2d<grid_index_type>(0, v.nX, 0, v.nY));
}

void GridBoundaryCondition_W_FreeSlip::enforce(Grid3<real,real>& w) const {
	Grid3_ZDirichlet<real, real>(w, 0)(blocked_range2d<grid_index_type>(0, w.nX, 0, w.nY));
	Grid3_XNeumann<real, real> wx(w);	wx(blocked_range2d<grid_index_type>(0, w.nY, 0, w.nZ));
	Grid3_YNeumann<real, real> wy(w);	wy(blocked_range2d<grid_index_type>(0, w.nX, 0, w.nZ));
}

// correct
void GridBoundaryCondition_V_CeilingSeparation::enforce(Grid2<real,real>& v) const {
	for ( Grid2r::index_type i = 0; i < v.nX; ++i )	{
		v(i,v.nY-1) = -ratio * v.dy;
		//v(i,nY-1) = -ratio * h;		
	}
}

void GridBoundaryCondition_V_CeilingSeparation::enforce(Grid3<real,real>& v) const {
	for ( Grid2r::index_type i = 0; i < v.nX; ++i )	
		for( Grid2r::index_type k = 0; k < v.nZ; ++k)	
			v(i,v.nY-1,k) = -ratio * v.dy;
			//v(i,nY,k) = -ratio * dy;			
}