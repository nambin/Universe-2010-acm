#pragma once

namespace fluidspace {

template<typename TRACER> 
void CFDAdvection_SemiLagrangian_SRL<TRACER>::post_process() const {
	log_info << std::endl << "### Log for CFDAdvection_SemiLagrangian_SRL _ Post Process ###" << std::endl;

	log_info << "Time for total spent time advection process       _ " << time_total << std::endl;
	log_info << "Time for total spent time only for pure advection _ " << time_only_advection << std::endl << std::endl;

#ifdef USE_OPENMP
	log_info << "CAUTION !!! : If OpenMP is Used _ Clipping Percentage information cannot be believed !!!!!!!!!!!!!!!!" << std::endl;
#endif
	log_info << "Number of Total Interpolation : " << num_itrp << std::endl;
	log_info << "Number of Clipping            : " << num_clipping << std::endl;
	log_info << "Percentage of Clipping        : " << float(num_clipping) / float(num_itrp) * 100.0f << std::endl;
}

// Need to Parallelize SRL
template<typename TRACER> template<typename CFD> 
void CFDAdvection_SemiLagrangian_SRL<TRACER>::advect_srl(InterfaceTracking_SRL2D& srl, const CFD& field, real dt) const {

	basis::Timer timer_all, timer_adv;

	AABB2D<real> aabb = srl.srl->getAABB();	

	timer_adv.reset();
	CFDAdvection_SemiLagrangian<TRACER>::advect(srl.srl0->phi, srl.srl->phi, field, dt);	
	time_only_advection += timer_adv.getTime();

	// Extract Narrow-band Cells using Coarse Grid Phi
	SRL2D<real>::extract_narrowband(srl.srl0->indices, srl.srl0->phi, srl.srl0->narrowbandWidth, unsigned int(SRL<real>::NEED_REFINE), unsigned int(SRL<real>::NOT_REFINE));

	// Add Narrow-band Cells which have sub-cell interfaces
	srl.add_narrowband_subcell_interface(srl.srl0->indices, srl.srl->nodes, unsigned int(SRL<real>::NEED_REFINE));

	// Subtract Narrow-band Cells which is contained in Kinematic Boundary
	srl.subtract_narrowband_kinematicbdry(srl.srl0->indices, unsigned int(SRL<real>::NOT_REFINE));

	// Refine Grid
	srl.srl0->refine_subgrid(srl.srl0->indices, unsigned int(SRL<real>::NEED_REFINE), unsigned int(SRL<real>::NOT_REFINE));

	SRL<real>::SRL_Interpolator interp(clipping);
	SRL<real>::SRL_Lerp lerp;
	SRL<real>::SRL_Coarse_Lerp coarse_lerp;

	timer_adv.reset();

#ifdef USE_OPENMP	
	int my_end = int(srl.srl0->nodes.size());
	int i;

#pragma omp parallel for
#else	
	std::size_t my_end = srl.srl0->nodes.size();
	std::size_t i;
#endif
	for(i = 0; i < my_end; ++i ) {		
		index_type ii, jj;
		Vector2<real> pos, posNew, x;
		Vector2<unsigned int> rs, re;

		SRLNode2D<real>& node = srl.srl0->nodes[i];
		srl.srl0->get_updating_subgrid_index(rs, re, node.index);

		Vector2<real> pos0 = srl.srl0->phi.pos(node.index);

		for(unsigned int p=rs[0];p<re[0];++p) {
			for(unsigned int q=rs[1];q<re[1];++q) {
				pos = pos0 + Vector2<real>(srl.srl0->pos_sub[p], srl.srl0->pos_sub[q]);
				posNew = TRACER::backtrack(pos, field, aabb, dt);
				aabb.constrainedBdry(posNew);

				srl.srl->get_idx_normalized_coord(ii, jj, x, posNew);

				if(srl.srl->is_exist_node(ii,jj)) {
					if(srl.srl->uniform_quadrature == false)	node.val[p][q] = interp(*srl.srl, ii, jj, x);
					else										node.val[p][q] = lerp(*srl.srl, ii, jj, x);
				}
				else {
					node.val[p][q] = coarse_lerp(*srl.srl, ii, jj, x);
				}
			}
		}		
	}

	time_only_advection += timer_adv.getTime();

	srl.srl0->synchronize_grid_values();
	swap(srl.srl, srl.srl0);

	time_total += timer_all.getTime();
	num_itrp += interp.num_itrp;
	num_clipping += interp.num_clipping;
}

// Need to Parallelize SRL
template<typename TRACER> template<typename CFD> 
void CFDAdvection_SemiLagrangian_SRL<TRACER>::advect_srl(InterfaceTracking_SRL3D& srl, const CFD& field, real dt) const {

	basis::Timer timer_all, timer_adv;

	AABB<real> aabb = srl.srl->getAABB();	

	timer_adv.reset();
	CFDAdvection_SemiLagrangian<TRACER>::advect(srl.srl0->phi, srl.srl->phi, field, dt);
	time_only_advection += timer_adv.getTime();

	// Extract Narrow-band Cells using Coarse Grid Phi
	SRL3D<real>::extract_narrowband(srl.srl0->indices, srl.srl0->phi, srl.srl0->narrowbandWidth, unsigned int(SRL<real>::NEED_REFINE), unsigned int(SRL<real>::NOT_REFINE));

	// Add Narrow-band Cells which have sub-cell interfaces
	srl.add_narrowband_subcell_interface(srl.srl0->indices, srl.srl->nodes, unsigned int(SRL<real>::NEED_REFINE));

	// Subtract Narrow-band Cells which is contained in Kinematic Boundary
	srl.subtract_narrowband_kinematicbdry(srl.srl0->indices, unsigned int(SRL<real>::NOT_REFINE));

	// Refine Grid
	srl.srl0->refine_subgrid(srl.srl0->indices, unsigned int(SRL<real>::NEED_REFINE), unsigned int(SRL<real>::NOT_REFINE));

	SRL<real>::SRL_Interpolator interp(clipping);
	SRL<real>::SRL_Lerp lerp;
	SRL<real>::SRL_Coarse_Lerp coarse_lerp;

	timer_adv.reset();

#ifdef USE_OPENMP	
	int my_end = int(srl.srl0->nodes.size());
	int i;

#pragma omp parallel for
#else	
	std::size_t my_end = srl.srl0->nodes.size();
	std::size_t i;
#endif
	for(i = 0; i < my_end; ++i ) {		
		index_type ii, jj, kk;
		Vector3<real> pos, posNew, x;
		Vector3<unsigned int> rs, re;

		SRLNode3D<real>& node = srl.srl0->nodes[i];
		srl.srl0->get_updating_subgrid_index(rs, re, node.index);

		Vector3<real> pos0 = srl.srl0->phi.pos(node.index);

		for(unsigned int p=rs[0];p<re[0];++p) {
			for(unsigned int q=rs[1];q<re[1];++q) {
				for(unsigned int r=rs[2];r<re[2];++r) {
					pos = pos0 + Vector3<real>(srl.srl0->pos_sub[p], srl.srl0->pos_sub[q], srl.srl0->pos_sub[r]);
					posNew = TRACER::backtrack(pos, field, aabb, dt);
					aabb.constrainedBdry(posNew);

					srl.srl->get_idx_normalized_coord(ii, jj, kk, x, posNew);

					if(srl.srl->is_exist_node(ii,jj,kk)) {
						if(srl.srl->uniform_quadrature == false)	
							node.val[p][q][r] = interp(*srl.srl, ii, jj, kk, x);
						else										
							node.val[p][q][r] = lerp(*srl.srl, ii, jj, kk, x);
					}
					else {
						node.val[p][q][r] = coarse_lerp(*srl.srl, ii, jj, kk, x);
					}

				}				
			}
		}		
	}

	time_only_advection += timer_adv.getTime();

	srl.srl0->synchronize_grid_values();
	swap(srl.srl, srl.srl0);

	time_total += timer_all.getTime();
	num_itrp += interp.num_itrp;
	num_clipping += interp.num_clipping;
}

};		// fluidspace