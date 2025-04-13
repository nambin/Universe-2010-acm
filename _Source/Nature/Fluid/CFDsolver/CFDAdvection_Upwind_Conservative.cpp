#include <stdafx.h>
#include <Nature/Fluid/CFDsolver/CFDAdvection.h>

using namespace std;
using namespace fluidspace;

void CFDAdvection_Upwind_Conservative::advect(Grid2r& s, const Grid2r& s0, const CFDsolver2D& cfd, real dt, MacGrid2D<real>::FIELD_TYPE type, const GridBdry& b) const {
	advect_ee(s,s0,cfd,dt,b);
}

void CFDAdvection_Upwind_Conservative::advect(Grid3r& s, const Grid3r& s0, const CFDsolver3D& cfd, real dt, MacGrid3D<real>::FIELD_TYPE type, const GridBdry& b) const {
	advect_ee(s,s0,cfd,dt,b);
}

void CFDAdvection_Upwind_Conservative::advect(Grid2r& s, const Grid2r& s0, const VelocityField2D<real>& vel, real dt, const GridBdry& b) const {
	advect_ee(s,s0,vel,dt,b);
}

void CFDAdvection_Upwind_Conservative::advect(Grid3r& s, const Grid3r& s0, const VelocityField3D<real>& vel, real dt, const GridBdry& b) const {
	advect_ee(s,s0,vel,dt,b);
}