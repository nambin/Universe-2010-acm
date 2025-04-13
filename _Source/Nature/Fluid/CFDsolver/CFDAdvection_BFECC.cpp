#include <stdafx.h>
#include <Nature/Fluid/CFDsolver/CFDAdvection.h>
#include <Principle/Geometry/MacGrid/MacGridFunctor.h>

using namespace std;
using namespace arithmetic;
using namespace geometryspace;
using namespace fluidspace;

void CFDAdvection_BFECC::advect(Grid2r& s, const Grid2r& s0, const CFDsolver2D& cfd, real dt, MacGrid2D<real>::FIELD_TYPE type, const GridBdry& b) const {
	Grid2r sTmp;
	sTmp.resize(s);

	CFDAdvection_SemiLagrangian<>::advect_sl( s, s0, cfd, dt, type, b );	
	CFDAdvection_SemiLagrangian<>::advect_sl( sTmp, s, cfd, -dt, type, b );	

#ifdef USE_OPENMP
	int my_end = int(s.size());
	int i;

#pragma omp parallel for
#else
	index_type my_end = s.size();
	index_type i;
#endif
	for( i = 0; i < my_end; ++i ) 
		s[i] = - real(0.5) * sTmp[i] + real(1.5) * s0[i];

	CFDAdvection_SemiLagrangian<>::advect_sl( sTmp, s, cfd, dt, type, b );

	if(limiting) {
		s = sTmp;
		CFDAdvection_SemiLagrangian<>::limiting_sl(s, s0, cfd, dt, type);

		b.enforce(s);
	}	
	else
		s = sTmp;
}

void CFDAdvection_BFECC::advect(Grid3r& s, const Grid3r& s0, const CFDsolver3D& cfd, real dt, MacGrid3D<real>::FIELD_TYPE type, const GridBdry& b) const {
	Grid3r sTmp;
	sTmp.resize(s);

	CFDAdvection_SemiLagrangian<>::advect_sl( s, s0, cfd, dt, type, b );	
	CFDAdvection_SemiLagrangian<>::advect_sl( sTmp, s, cfd, -dt, type, b );	

#ifdef USE_OPENMP
	int my_end = int(s.size());
	int i;

#pragma omp parallel for
#else
	index_type my_end = s.size();
	index_type i;
#endif
	for( i = 0; i < my_end; ++i ) 
		s[i] = - real(0.5) * sTmp[i] + real(1.5) * s0[i];

	CFDAdvection_SemiLagrangian<>::advect_sl( sTmp, s, cfd, dt, type, b );

	if(limiting) {
		s = sTmp;
		CFDAdvection_SemiLagrangian<>::limiting_sl(s, s0, cfd, dt, type);

		b.enforce(s);
	}	
	else
		s = sTmp;
}

void CFDAdvection_BFECC::advect(Grid2r& s, const Grid2r& s0, const VelocityField2D<real>& vel, real dt, const GridBdry& b) const {
	Grid2r sTmp;
	sTmp.resize(s);
	
	CFDAdvection_SemiLagrangian<>::advect_sl( s, s0, vel, dt, b );	
	CFDAdvection_SemiLagrangian<>::advect_sl( sTmp, s, vel, -dt, b );	

#ifdef USE_OPENMP
	int my_end = int(s.size());
	int i;

#pragma omp parallel for
#else
	index_type my_end = s.size();
	index_type i;
#endif
	for( i = 0; i < my_end; ++i ) 
		s[i] = - real(0.5) * sTmp[i] + real(1.5) * s0[i];

	CFDAdvection_SemiLagrangian<>::advect_sl( sTmp, s, vel, dt, b );

	if(limiting) {
		s = sTmp;
		CFDAdvection_SemiLagrangian<>::limiting_sl(s, s0, vel, dt);

		b.enforce(s);
	}	
	else
		s = sTmp;
}	

void CFDAdvection_BFECC::advect(Grid3r& s, const Grid3r& s0, const VelocityField3D<real>& vel, real dt, const GridBdry& b) const {
	Grid3r sTmp;
	sTmp.resize(s);

	CFDAdvection_SemiLagrangian<>::advect_sl( s, s0, vel, dt, b );	
	CFDAdvection_SemiLagrangian<>::advect_sl( sTmp, s, vel, -dt, b );	

#ifdef USE_OPENMP
	int my_end = int(s.size());
	int i;

#pragma omp parallel for
#else
	index_type my_end = s.size();
	index_type i;
#endif
	for( i = 0; i < my_end; ++i ) 
		s[i] = - real(0.5) * sTmp[i] + real(1.5) * s0[i];

	CFDAdvection_SemiLagrangian<>::advect_sl( sTmp, s, vel, dt, b );

	if(limiting) {
		s = sTmp;
		CFDAdvection_SemiLagrangian<>::limiting_sl(s, s0, vel, dt);

		b.enforce(s);
	}	
	else
		s = sTmp;
}
