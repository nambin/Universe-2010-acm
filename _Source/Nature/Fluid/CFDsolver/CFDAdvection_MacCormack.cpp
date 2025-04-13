#include <stdafx.h>
#include <Nature/Fluid/CFDsolver/CFDAdvection.h>

using namespace std;
using namespace fluidspace;

void CFDAdvection_MacCormack::advect(Grid2r& s, const Grid2r& s0, const CFDsolver2D& cfd, real dt, MacGrid2D<real>::FIELD_TYPE type, const GridBdry& b) const {
	Grid2r sTmp;
	sTmp.resize(s);

	CFDAdvection_SemiLagrangian<>::advect_sl( s, s0, cfd, dt, type, b );	
	CFDAdvection_SemiLagrangian<>::advect_sl( sTmp, s, cfd, -dt, type, b );	

	for(index_type i = 0; i < s.size(); ++i)
		sTmp[i] = s[i] + 0.5f * (s0[i] - sTmp[i]);
	
	if(limiting) {
		if(revertsl)	CFDAdvection_SemiLagrangian<>::revert_sl(s, s0, sTmp, cfd, dt, type);
		else	{		s = sTmp; CFDAdvection_SemiLagrangian<>::limiting_sl(s, s0, cfd, dt, type); }
	}	
	else
		s = sTmp;
}

void CFDAdvection_MacCormack::advect(Grid3r& s, const Grid3r& s0, const CFDsolver3D& cfd, real dt, MacGrid3D<real>::FIELD_TYPE type, const GridBdry& b) const {

	Grid3r sTmp;
	sTmp.resize(s);

	CFDAdvection_SemiLagrangian<>::advect_sl( s, s0, cfd, dt, type, b );	
	CFDAdvection_SemiLagrangian<>::advect_sl( sTmp, s, cfd, -dt, type, b );	

	for(index_type i = 0; i < s.size(); ++i)
		sTmp[i] = s[i] + 0.5f * (s0[i] - sTmp[i]);

	if(limiting) {
		if(revertsl)	CFDAdvection_SemiLagrangian<>::revert_sl(s, s0, sTmp, cfd, dt, type);
		else	{		s = sTmp; CFDAdvection_SemiLagrangian<>::limiting_sl(s, s0, cfd, dt, type); }
	}	
	else
		s = sTmp;
}

void CFDAdvection_MacCormack::advect(Grid2r& s, const Grid2r& s0, const VelocityField2D<real>& vel, real dt, const GridBdry& b) const {
	Grid2r sTmp;
	sTmp.resize(s);

	CFDAdvection_SemiLagrangian<>::advect_sl( s, s0, vel, dt, b );	
	CFDAdvection_SemiLagrangian<>::advect_sl( sTmp, s, vel, -dt, b );		

	for(index_type i = 0; i < s.size(); ++i)
		sTmp[i] = s[i] + 0.5f * (s0[i] - sTmp[i]);

	if(limiting) {
		if(revertsl)	CFDAdvection_SemiLagrangian<>::revert_sl(s, s0, sTmp, vel, dt);
		else	{		s = sTmp; CFDAdvection_SemiLagrangian<>::limiting_sl(s, s0, vel, dt); }
	}	
	else
		s = sTmp;
}	

void CFDAdvection_MacCormack::advect(Grid3r& s, const Grid3r& s0, const VelocityField3D<real>& vel, real dt, const GridBdry& b) const {
	Grid3r sTmp;
	sTmp.resize(s);

	CFDAdvection_SemiLagrangian<>::advect_sl( s, s0, vel, dt, b );	
	CFDAdvection_SemiLagrangian<>::advect_sl( sTmp, s, vel, -dt, b );		

	for(index_type i = 0; i < s.size(); ++i)
		sTmp[i] = s[i] + 0.5f * (s0[i] - sTmp[i]);

	if(limiting) {
		if(revertsl)	CFDAdvection_SemiLagrangian<>::revert_sl(s, s0, sTmp, vel, dt);
		else	{		s = sTmp; CFDAdvection_SemiLagrangian<>::limiting_sl(s, s0, vel, dt); }
	}	
	else
		s = sTmp;
}