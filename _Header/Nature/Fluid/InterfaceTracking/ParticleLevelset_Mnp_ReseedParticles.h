
template<typename PLS, typename HEAPGRID, typename IDXVEC>
void PLS_ParticleMnp::collectParticles(HEAPGRID& heaps, PLS& pls) const {

	typedef typename HEAPGRID::value_type		particleHeap;
	typedef typename particleHeap::value_type	particlePriority;
	
	IDXVEC idx;
	for(size_t i=0;i<pls.particles.size();++i) {

		if(pls.particles[i].isDeleted)
			continue;

		real _phi_p = pls.phi(pls.particles[i].pos);

		if(abs(_phi_p) > b_max && !pls.particles[i].isEscaped) {
			pls.particles[i].isDeleted = true;
			continue;
		}

		pls.phi.getClosestIdx(idx, pls.particles[i].pos);

		if(!pls.particles[i].isEscaped) {
			particlePriority pp(&pls.particles[i], pls.particles[i].sign * _phi_p - pls.particles[i].radius);
			particleHeap &heap = heaps(idx);

			if(heap.size() < numParticlesPerCell) 
				heap.push(pp);
			else {
				const particlePriority &top_pp = heap.top();
				if(top_pp.second > pp.second) {
					top_pp.first->isDeleted = true;
					heap.pop();
					heap.push(pp);
				}
				else
					pls.particles[i].isDeleted = true;
			}
		}
	}	
}