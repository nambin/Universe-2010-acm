#include <stdafx.h>
#include <Nature/Fluid/InterfaceTracking/ParticleLevelset_ParticleAdvect.h>
#include <Nature/Fluid/CFDsolver/CFDTracking.h>

using namespace std;
using namespace fluidspace;
	
void PLS_PartlcleAdvector_RK2::advect(InterfaceTracking_PLS2D& pls, const CFDsolver2D& cfd, real dt) const {

	real Lx = cfd.Lx;
	real Ly = cfd.Ly;
	real eps = InterfaceTracking_PLS2D::Epsilon();
	
	// 2nd RK
	for (size_t i=0;i<pls.particles.size();++i) {
		if(pls.particles[i].isDeleted)
			continue;

		pls.particles[i].pos = CFDTracking_RK2::forwardtrack(pls.particles[i].pos, cfd, dt);

		// Boundary clipping
		if ( pls.particles[i].pos[0] <= 0.0f )		pls.particles[i].pos[0] = eps;
		else if ( pls.particles[i].pos[0] >= Lx )	pls.particles[i].pos[0] = Lx - eps;
		if ( pls.particles[i].pos[1] <= 0.0f )		pls.particles[i].pos[1] = eps;
		else if ( pls.particles[i].pos[1] >= Ly )	pls.particles[i].pos[1] = Ly - eps;
	}
}

void PLS_PartlcleAdvector_RK2::advect(InterfaceTracking_PLS3D& pls, const CFDsolver3D& cfd, real dt) const {

	real Lx = cfd.Lx;
	real Ly = cfd.Ly;
	real Lz = cfd.Lz;
	real eps = InterfaceTracking_PLS3D::Epsilon();

//#ifdef USE_OPENMP	
//	int my_end = int(pls.particles.size());
//	int i;
//
//#pragma omp parallel for
//#else	
	std::size_t my_end = pls.particles.size();
	std::size_t i;
//#endif
	for( i = 0; i < my_end; ++i ) {		
		if(pls.particles[i].isDeleted)
			continue;

		pls.particles[i].pos = CFDTracking_RK2::forwardtrack(pls.particles[i].pos, cfd, dt);

		// Boundary clipping
		if ( pls.particles[i].pos[0] <= 0.0f )		pls.particles[i].pos[0] = eps;
		else if ( pls.particles[i].pos[0] >= Lx )	pls.particles[i].pos[0] = Lx - eps;
		if ( pls.particles[i].pos[1] <= 0.0f )		pls.particles[i].pos[1] = eps;
		else if ( pls.particles[i].pos[1] >= Ly )	pls.particles[i].pos[1] = Ly - eps;
		if ( pls.particles[i].pos[2] <= 0.0f )		pls.particles[i].pos[2] = eps;
		else if ( pls.particles[i].pos[2] >= Lz )	pls.particles[i].pos[2] = Lz - eps;
	}
}

void PLS_PartlcleAdvector_RK2::advect(InterfaceTracking_PLS2D& pls, const VelocityField2D<real>& vel, real dt) const {
	real Lx = pls.Lx;
	real Ly = pls.Ly;
	real eps = InterfaceTracking_PLS2D::Epsilon();

	// 2nd RK
	for (size_t i=0;i<pls.particles.size();++i) {
		if(pls.particles[i].isDeleted)
			continue;

		Vector2<real> mid = pls.particles[i].pos + (0.5f * dt) * vel.getVel(pls.particles[i].pos);
		pls.constrainedBdry(mid);
		
		pls.particles[i].pos += (dt * vel.getVel(mid));

		// Boundary clipping
		if ( pls.particles[i].pos[0] <= 0.0f )		pls.particles[i].pos[0] = eps;
		else if ( pls.particles[i].pos[0] >= Lx )	pls.particles[i].pos[0] = Lx - eps;
		if ( pls.particles[i].pos[1] <= 0.0f )		pls.particles[i].pos[1] = eps;
		else if ( pls.particles[i].pos[1] >= Ly )	pls.particles[i].pos[1] = Ly - eps;
	}
}

void PLS_PartlcleAdvector_RK2::advect(InterfaceTracking_PLS3D& pls, const VelocityField3D<real>& vel, real dt) const {
	real Lx = pls.Lx;
	real Ly = pls.Ly;
	real Lz = pls.Lz;
	real eps = InterfaceTracking_PLS3D::Epsilon();

#ifdef USE_OPENMP	
	int my_end = int(pls.particles.size());
	int i;

#pragma omp parallel for
#else	
	std::size_t my_end = pls.particles.size();
	std::size_t i;
#endif
	for( i = 0; i < my_end; ++i ) {	
		if(pls.particles[i].isDeleted)
			continue;

		Vector3<real> mid = pls.particles[i].pos + (0.5f * dt) * vel.getVel(pls.particles[i].pos);
		pls.constrainedBdry(mid);

		pls.particles[i].pos += (dt * vel.getVel(mid));

		// Boundary clipping
		if ( pls.particles[i].pos[0] <= 0.0f )		pls.particles[i].pos[0] = eps;
		else if ( pls.particles[i].pos[0] >= Lx )	pls.particles[i].pos[0] = Lx - eps;
		if ( pls.particles[i].pos[1] <= 0.0f )		pls.particles[i].pos[1] = eps;
		else if ( pls.particles[i].pos[1] >= Ly )	pls.particles[i].pos[1] = Ly - eps;
		if ( pls.particles[i].pos[2] <= 0.0f )		pls.particles[i].pos[2] = eps;
		else if ( pls.particles[i].pos[2] >= Lz )	pls.particles[i].pos[2] = Lz - eps;
	}
}