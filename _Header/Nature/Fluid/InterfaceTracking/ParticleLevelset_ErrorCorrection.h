#ifndef _FLUIDSPACE_PARTICLE_LEVELSET_ERROR_CORRECTION_H_
#define _FLUIDSPACE_PARTICLE_LEVELSET_ERROR_CORRECTION_H_

#include <vector>
#include "../Environment/Random.h"
#include "../InterfaceTracking/InterfaceTracking_PLS2D.h"
#include "../InterfaceTracking/InterfaceTracking_PLS3D.h"

namespace fluidspace {

class NATUREDLL PLS_ErrorCorrection {
public : 
	void errorCorrect(InterfaceTracking_PLS2D&, bool setNeumann = true) const;	
	void errorCorrect(InterfaceTracking_PLS3D&, bool setNeumann = true) const;	

protected : 
	template<typename PLS, typename PHIGRID, typename BOOLGRID, typename IDXVEC>
	static void getCorrectedPhi(PHIGRID& pos, PHIGRID& neg, BOOLGRID& isModified, const PLS&, const std::vector<IDXVEC>&);		// !!!
};

template<typename PLS, typename PHIGRID, typename BOOLGRID, typename IDXVEC>
void PLS_ErrorCorrection::getCorrectedPhi(PHIGRID& phiPos, PHIGRID& phiNeg, BOOLGRID& isModified, const PLS& pls, const std::vector<IDXVEC>& nIdx) {
	
	IDXVEC idx, _idx;
	for(size_t idxE=0;idxE<pls.idxEscapedParticles.size();++idxE) {
		typename const PLS::Particle& particle = pls.particles[pls.idxEscapedParticles[idxE]];

//	for(size_t i=0;i<pls.particles.size();++i) {
//		typename const PLS::Particle& particle = pls.particles[i];

//		real _phi = pls.phi(particle.pos);
//		if(_phi * particle.sign > 0 && abs(_phi) > particle.radius)
		//if(_phi * particle.sign > 0 || abs(_phi) < particle.radius)
//			continue;

		pls.phi.getIdx(idx,particle.pos);
		if(particle.sign == 1) {
			for(size_t i=0;i<nIdx.size();++i) {
				_idx = idx + nIdx[i];
				phiPos(_idx) = MAX(phiPos(_idx), particle.getPhi(pls.phi.pos(_idx)));
				isModified(_idx) = true;
			}			
		}
		else {
			for(size_t i=0;i<nIdx.size();++i) {
				_idx = idx + nIdx[i];
				phiNeg(_idx) = MIN(phiNeg(_idx), particle.getPhi(pls.phi.pos(_idx)));
				isModified(_idx) = true;
			}
		}	
	}
}

};		// fluidspace

#endif