
template<typename PLS, typename NORMALGRID, typename IDXVEC>	
void PLS_ParticleMnp::add_particles(PLS& pls, const IDXVEC& idx, size_t num, const NORMALGRID& normals) const {

	typedef typename NORMALGRID::value_type POSVEC; 

	if(abs(pls.phi(idx)) > b_max)
		return;

	for(size_t k=0;k<num;++k) {		
		POSVEC posParticle = pls.phi.pos(idx) + (POSVEC(random) - POSVEC(0.5f)) * pls.h;
		typename PLS::Particle particle(posParticle, pls.phi(posParticle));

		real phigoal = (b_max - b_min) * random() + b_min;
		phigoal *= particle.sign;

		attractParticle(particle, phigoal, normals, pls);
		if(particle.isDeleted) {
			--k;
			continue;
		}

		pls.particles.push_back(particle);
	}
}

template<typename PLS>
void PLS_ParticleMnp::removeDeepParticles(PLS& pls) const {
	real _phi;
	for(size_t i=0;i<pls.particles.size();++i) {
		if(pls.particles[i].isDeleted)
			continue;

		_phi = pls.phi(pls.particles[i].pos);

		if(SIGN(_phi) == pls.particles[i].sign && abs(_phi) > b_max)
			pls.particles[i].isDeleted = true;
	}
}

template<typename PLS>
void PLS_ParticleMnp::delete_isDeletedParticles(PLS& pls) const {
	vector<typename PLS::Particle> _tmpParticles = pls.particles;
	pls.particles.clear();
	for(size_t i=0;i<_tmpParticles.size();++i) 
		if(!_tmpParticles[i].isDeleted)
			pls.particles.push_back(_tmpParticles[i]);
}
