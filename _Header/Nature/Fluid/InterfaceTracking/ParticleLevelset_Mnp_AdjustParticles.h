
template<typename PLS, typename NORMALGRID>	
void PLS_ParticleMnp::attractParticle(typename PLS::Particle& particle, real phigoal, const NORMALGRID &normals, const PLS& pls) const {

	typedef typename NORMALGRID::value_type POSVEC; 

	const size_t nIter = 15;

	real _phi, _phiNewAbs, lambda;
	POSVEC _normal, newPos;

	for(size_t iter=0;iter<nIter;++iter) {

		lambda = 1.0f;
		_phi = pls.phi(particle.pos);
		_normal = normals(particle.pos);

		newPos = particle.pos + lambda * (phigoal - _phi) * _normal;
		while(true) {
			if(!pls.isInsideDomain(newPos)) {
				lambda *= 0.5f;
				newPos = particle.pos + lambda * (phigoal - _phi) * _normal;
			}
			else
				break;
		}

		_phiNewAbs = abs(pls.phi(newPos));
		if(lambda != 1.0f && (_phiNewAbs < b_min || _phiNewAbs > b_max)) {
			lambda *= 0.5f;
			newPos = particle.pos + lambda * (phigoal - _phi) * _normal;
		}
		else if(lambda != 1.0f) {
			particle.pos = newPos;
			break;				
		}

		particle.pos = newPos;
	}

	_phi = pls.phi(particle.pos);
	if(abs(_phi) < b_min || abs(_phi) > b_max) 
		particle.isDeleted = true;
	else {
		particle.sign = SIGN(_phi);
		particle.radius = MIN(MAX(abs(_phi),r_min),r_max);
	}	
}

template<typename PLS>
void PLS_ParticleMnp::adjustRadius(PLS& pls) const {
	detect_EscapedParticles(pls);

	for(size_t i=0;i<pls.particles.size();++i) {
		if(pls.particles[i].isDeleted)
			continue;

		if(pls.particles[i].isEscaped) {
			pls.particles[i].radius = r_min;		// [Enright et al. 2004]
		}
		else {
			const real _phi = pls.phi(pls.particles[i].pos);

			if(SIGN(_phi) != pls.particles[i].sign)
				pls.particles[i].radius = r_min;
			else
				pls.particles[i].radius = MIN(MAX(abs(_phi),r_min),r_max);
		}
	}
}
