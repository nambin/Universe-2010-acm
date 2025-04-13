#ifndef _FLUIDSPACE_PARTICLE_LEVELSET_MANIPULATOR_H_
#define _FLUIDSPACE_PARTICLE_LEVELSET_MANIPULATOR_H_

#include <queue>
#include "../Environment/Random.h"
#include "../InterfaceTracking/InterfaceTracking_PLS2D.h"
#include "../InterfaceTracking/InterfaceTracking_PLS3D.h"

namespace fluidspace {

class NATUREDLL PLS_ParticleMnp {
public : 
	PLS_ParticleMnp();

	void initialize(InterfaceTracking_PLS2D&) const;
	void initialize(InterfaceTracking_PLS3D&) const;

	void reseed(InterfaceTracking_PLS2D&) const;
	void reseed(InterfaceTracking_PLS3D&) const;

	template<typename PLS>
	void removeDeepParticles(PLS&) const;

	template<typename PLS>
	void adjustRadius(PLS&) const;																					// !!!

	template<typename PLS>
	static void detect_EscapedParticles(PLS&);

protected : 	
	template<typename PLS, typename NORMALGRID, typename IDXVEC>	
	void add_particles(PLS&, const IDXVEC&, size_t num, const NORMALGRID& normals) const;

	template<typename PLS, typename NORMALGRID>	
	void attractParticle(typename PLS::Particle&, real phigoal, const NORMALGRID &normals, const PLS&) const;		// !!!

	template<typename PLS, typename HEAPGRID, typename IDXVEC>
	void collectParticles(HEAPGRID&, PLS&) const;																	// !!!

	template<typename PLS, typename HEAPGRID, typename NORMALGRID, typename IDXVEC>
	void reseeding(PLS&, const HEAPGRID&) const;

	template<typename PLS>
	void delete_isDeletedParticles(PLS&) const;

public : 
	real				b_min_frac, b_max_frac;
	real				r_min_frac, r_max_frac;
	size_t				numParticlesPerCell_1D;
	size_t				reseed_cycle;
	
protected : 
	mutable real		b_min, b_max;
	mutable real		r_min, r_max;
	mutable size_t		numParticlesPerCell;	
};

#include "ParticleLevelset_Mnp_AddRemoveParticles.h"
#include "ParticleLevelset_Mnp_ReseedParticles.h"
#include "ParticleLevelset_Mnp_AdjustParticles.h"

template<typename PLS>
void PLS_ParticleMnp::detect_EscapedParticles(PLS& pls) {
	pls.idxEscapedParticles.clear();

	for(size_t i = 0;i<pls.particles.size();++i) {

		if(pls.particles[i].isDeleted)
			continue;

		const real _phi = pls.phi(pls.particles[i].pos);

		// test
		//typename PLS::Particle par = pls.particles[i];

		if(SIGN(_phi) != pls.particles[i].sign && abs(_phi) > pls.particles[i].radius) {
			pls.particles[i].isEscaped = true;
			pls.idxEscapedParticles.push_back(unsigned int(i));
		}
		else 
			pls.particles[i].isEscaped = false;
	}
} 

};			// fluidspace

#endif
