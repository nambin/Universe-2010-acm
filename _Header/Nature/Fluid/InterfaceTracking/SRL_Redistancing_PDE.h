#pragma once

namespace fluidspace {

template<typename DEV_SRL, typename DEV_LV> 
void SRL_Redistancing_PDE::redistance_pde(InterfaceTracking_SRL2D& srl) const {

	Timer timer_all, timer_fine;

	typedef real T;

	SRL2D<real>* src = srl.srl, *dst = srl.srl0;

	// Re-distancing Coarse Grid
	src->modify_coarse_phi_except_near_interface();
	geometryspace::Levelset_Redistancing_PDE<real>::redistance_pde(dst->phi, src->phi, Redistancing_EE<Upwind_1st<Grid2<T,T>>>(*this));

	// Extract Narrow-band Cells using Coarse Grid Phi
	SRL2D<T>::extract_narrowband(dst->indices, dst->phi, src->narrowbandWidth, unsigned int(SRL<T>::NEED_REFINE), unsigned int(SRL<T>::NOT_REFINE));

	// Add Narrow-band Cells which have sub-cell interfaces
	srl.add_narrowband_subcell_interface(dst->indices, src->nodes, unsigned int(SRL<T>::NEED_REFINE));

	// Subtract Narrow-band Cells which is contained in Kinematic Boundary
	srl.subtract_narrowband_kinematicbdry(srl.srl0->indices, unsigned int(SRL<real>::NOT_REFINE));

	// Refine Grid
	dst->refine_subgrid(dst->indices, unsigned int(SRL<T>::NEED_REFINE), unsigned int(SRL<T>::NOT_REFINE));

	// Copy Values & Structures
	dst->set_value(*src);	dst->synchronize_grid_values();
	src->set(*dst);

//	if(narrowband_calculation_only) redistance_pde_narrowband<DEV_SRL, DEV_LV>(srl); 
//	else							redistance_pde_wholedomain<DEV_SRL, DEV_LV>(srl); 

	timer_fine.reset();
	redistance_pde_narrowband<DEV_SRL, DEV_LV>(srl); 
	time_only_fine_grid += timer_fine.getTime(); 

	time_total += timer_all.getTime();

	srl.srl->print_debug_info();
}

template<typename DEV_SRL, typename DEV_LV> 
void SRL_Redistancing_PDE::redistance_pde(InterfaceTracking_SRL3D& srl) const {
	typedef real T;

	Timer timer_all, timer_fine;
	Timer timer;

	SRL3D<real>* src = srl.srl, *dst = srl.srl0;

	// Re-distancing Coarse Grid
	src->modify_coarse_phi_except_near_interface();
	geometryspace::Levelset_Redistancing_PDE<real>::redistance_pde(dst->phi, src->phi, Redistancing_EE<Upwind_1st<Grid3<T,T>>>(*this));

	timer.logTime("1. Redist Coarse", log_info); timer.reset();

	// Extract Narrow-band Cells using Coarse Grid Phi
	SRL3D<T>::extract_narrowband(dst->indices, dst->phi, src->narrowbandWidth, unsigned int(SRL<T>::NEED_REFINE), unsigned int(SRL<T>::NOT_REFINE));

	timer.logTime("2. Extract Narrowband", log_info); timer.reset();

	// Add Narrow-band Cells which have sub-cell interfaces
	srl.add_narrowband_subcell_interface(dst->indices, src->nodes, unsigned int(SRL<T>::NEED_REFINE));

	timer.logTime("3. Add Narrowband", log_info); timer.reset();

	// Subtract Narrow-band Cells which is contained in Kinematic Boundary
	srl.subtract_narrowband_kinematicbdry(srl.srl0->indices, unsigned int(SRL<real>::NOT_REFINE));

	timer.logTime("4. Subtract Narrowband inside kinematic bdry", log_info); timer.reset();

	// Refine Grid
	dst->refine_subgrid(dst->indices, unsigned int(SRL<T>::NEED_REFINE), unsigned int(SRL<T>::NOT_REFINE));

	timer.logTime("5. Refine Subgrid", log_info); timer.reset();

	// Copy Values & Structures
	dst->set_value(*src);	dst->synchronize_grid_values();
	src->set(*dst);

	timer.logTime("6. Copy Grid Values Subgrid", log_info); timer.reset();

//	if(narrowband_calculation_only) redistance_pde_narrowband<DEV_SRL, DEV_LV>(srl); 
//	else							redistance_pde_wholedomain<DEV_SRL, DEV_LV>(srl); 

	timer_fine.reset();
	redistance_pde_narrowband<DEV_SRL, DEV_LV>(srl);
	time_only_fine_grid += timer_fine.getTime(); 

	time_total += timer_all.getTime();

	timer.logTime("7. Redist Fine Grid", log_info); timer.reset();

	srl.srl->print_debug_info();
}

template<typename DEV_SRL, typename DEV_LV> 
void SRL_Redistancing_PDE::redistance_pde_wholedomain(InterfaceTracking_SRL2D& srl) const {
	typedef real T;

	SRL2D<real>* src = srl.srl, *dst = srl.srl0;

	// Calculate time-step Tau & Required Iteration Number
	const T h = src->h_sub[0], h2 = h*h;
	const T tau = cfl * h;
	size_t nIter = MAX(size_t((src->narrowbandWidth + h) / tau) + 1, min_iter);		

	// Iteration
	for(size_t iter = 0; iter < nIter; ++iter) {

#ifdef USE_OPENMP
		int my_end = int(src->phi.nX-1);
		int i;

#pragma omp parallel for
#else
		index_type my_end = src->phi.nX-1;
		index_type i;
#endif

		for( i = 0; i < my_end ; ++i ) {		

#define CLIP(x,xp,xm) if((x)>0) { xp = (x); xm = 0; } else { xp = 0; xm  = (x); }

			DEV_SRL upwind(*src);
			DEV_LV upwind_coarse(src->phi);

			T a,b,c,d,sign;
			T ap,bp,cp,dp;
			T am,bm,cm,dm;

			Vector2<unsigned int> rs, re;

			index_type		coarse_idx[2];
			unsigned int	sub_idx[2]; 

			coarse_idx[0] = i;
			for( index_type j = 0; j < src->phi.nY-1 ; ++j ) {

				coarse_idx[1] = j;
				if(src->indices(i,j) == SRL<T>::NOT_REFINE) 
				{
					if(i == 0 || j == 0)
						continue;					

					a = upwind_coarse.dx(i,j,true);
					b = upwind_coarse.dx(i,j,false);
					c = upwind_coarse.dy(i,j,true);
					d = upwind_coarse.dy(i,j,false);

					CLIP(a,ap,am);	CLIP(b,bp,bm);
					CLIP(c,cp,cm);	CLIP(d,dp,dm);

					sign = sign_sussman94(src->phi(i,j), h2);

					if(sign > 0)	{ dst->phi(i,j) = src->phi(i,j) - tau * sign * (sqrt(ap*ap + bm*bm + cp*cp + dm*dm) - 1);	}
					else			{ dst->phi(i,j) = src->phi(i,j)	- tau * sign * (sqrt(am*am + bp*bp + cm*cm + dp*dp) - 1);	}

					continue;
				}

				const SRLNode2D<T>&	src_node = src->node(i,j);
				SRLNode2D<T>&		dst_node = dst->node(i,j);			

				src->get_updating_subgrid_index(rs, re, src_node.index);

				for(unsigned int p=rs[0];p<re[0];++p) {
					sub_idx[0] = p;
					for(unsigned int q=rs[1];q<re[1];++q) {
						sub_idx[1] = q;

						upwind.dx(a,b,src_node,sub_idx);
						upwind.dy(c,d,src_node,sub_idx);

						CLIP(a,ap,am);	CLIP(b,bp,bm);
						CLIP(c,cp,cm);	CLIP(d,dp,dm);

						sign = sign_sussman94(src_node(p,q), h2);
						//sign = sign_peng99(src_node(p,q), h2);		// Is it possible ?

						if(sign > 0)	{ dst_node(p,q) = src_node(p,q) - tau * sign * (sqrt(ap*ap + bm*bm + cp*cp + dm*dm) - 1); }
						else			{ dst_node(p,q) = src_node(p,q) - tau * sign * (sqrt(am*am + bp*bp + cm*cm + dp*dp) - 1); }
					}
				}
			}

#undef CLIP
		}

		swap(src, dst);		
		src->synchronize_grid_values();

	}		// iter

	srl.srl = src; 
	srl.srl0 = dst;
}

template<typename DEV_SRL, typename DEV_LV> 
void SRL_Redistancing_PDE::redistance_pde_narrowband(InterfaceTracking_SRL2D& srl) const {
	typedef real T;

	SRL2D<real>* src = srl.srl, *dst = srl.srl0;

	// Calculate time-step Tau & Required Iteration Number	
	const T h = src->h_sub[0], h2 = h*h;
	const T tau = cfl * h;
	size_t nIter = MAX(size_t((src->narrowbandWidth + h) / tau) + 1, min_iter);		

	const T eps = 1.0e-3 * h;

	// Iteration
	for(size_t iter = 0; iter < nIter; ++iter) {

#ifdef USE_OPENMP
		int my_end = int(src->nodes.size());
		int i;

#pragma omp parallel for
#else
		std::size_t my_end = src->nodes.size();
		std::size_t i;
#endif
		for( i = 0; i < my_end ; ++i ) {		

#define CLIP(x,xp,xm) if((x)>0) { xp = (x); xm = 0; } else { xp = 0; xm  = (x); }

			DEV_SRL upwind(*src);

			const SRLNode2D<T>&	src_node = src->nodes[i];
			SRLNode2D<T>&		dst_node = dst->nodes[i];	

			T a,b,c,d,sign;
			T ap,bp,cp,dp;
			T am,bm,cm,dm;

			Vector2<unsigned int> rs, re;

			index_type		coarse_idx[2];
			unsigned int	sub_idx[2]; 								

			coarse_idx[0] = src_node.index[0];
			coarse_idx[1] = src_node.index[1];	

			src->get_updating_subgrid_index(rs, re, src_node.index);

			for(unsigned int p=rs[0];p<re[0];++p) {
				sub_idx[0] = p;
				for(unsigned int q=rs[1];q<re[1];++q) {
					sub_idx[1] = q;

					upwind.dx(a,b,src_node,sub_idx);
					upwind.dy(c,d,src_node,sub_idx);

					CLIP(a,ap,am);	CLIP(b,bp,bm);
					CLIP(c,cp,cm);	CLIP(d,dp,dm);

//					sign = sign_sussman94(src_node(p,q), h2);
//					sign = sign_peng99(src_node(p,q), h2);		// Is it possible ?

					const T _phi = (src_node(p,q) == 0) ? - eps : src_node(p,q);
					const T denom2 = SQ(_phi) + h2 * ( SQ( ((abs(a) > abs(b)) ? a : b) ) 
													 + SQ( ((abs(c) > abs(d)) ? c : d) ) );
					sign =  T(_phi / sqrt(denom2));					

					if(sign > 0)	{ dst_node(p,q) = src_node(p,q) - tau * sign * (sqrt(ap*ap + bm*bm + cp*cp + dm*dm) - 1); if(dst_node(p,q) < 0) dst_node(p,q) = + eps; }
					else			{ dst_node(p,q) = src_node(p,q) - tau * sign * (sqrt(am*am + bp*bp + cm*cm + dp*dp) - 1); if(dst_node(p,q) > 0) dst_node(p,q) = - eps; }
				}
			}

#undef CLIP
		}

		swap(src, dst);		
		src->synchronize_grid_values();	

	}		// iter

	srl.srl = src; 
	srl.srl0 = dst;
}

template<typename DEV_SRL, typename DEV_LV> 
void SRL_Redistancing_PDE::redistance_pde_narrowband(InterfaceTracking_SRL3D& srl) const {
	typedef real T;

	SRL3D<real>* src = srl.srl, *dst = srl.srl0;

	// Calculate time-step Tau & Required Iteration Number
	const T h = src->h_sub[0], h2 = h*h;
	const T tau = cfl * h;
	size_t nIter = MAX(size_t((src->narrowbandWidth) / tau) + 1, min_iter);		

	const T eps = 1.0e-3 * h;

	log_info << "Number of Iteration & Narrowband for redistancing fine : " << nIter << " " << src->narrowband << std::endl; 
	log_info << "Number of Nodes : " << src->nodes.size() << std::endl;

	// Iteration
	for(size_t iter = 0; iter < nIter; ++iter) {

#ifdef USE_OPENMP
		int my_end = int(src->nodes.size());
		int i;

#pragma omp parallel for
#else
		std::size_t my_end = src->nodes.size();
		std::size_t i;
#endif
		for( i = 0; i < my_end ; ++i ) {		

#define CLIP(x,xp,xm) if((x)>0) { xp = (x); xm = 0; } else { xp = 0; xm  = (x); }

			DEV_SRL upwind(*src);

			const SRLNode3D<T>&	src_node = src->nodes[i];
			SRLNode3D<T>&		dst_node = dst->nodes[i];	

			T a,b,c,d,e,f,sign;
			T ap,bp,cp,dp,ep,fp;
			T am,bm,cm,dm,em,fm;

			Vector3<unsigned int> rs, re;

			index_type		coarse_idx[3];
			unsigned int	sub_idx[3]; 								

			coarse_idx[0] = src_node.index[0];
			coarse_idx[1] = src_node.index[1];	
			coarse_idx[2] = src_node.index[2];	

			src->get_updating_subgrid_index(rs, re, src_node.index);

			for(unsigned int p=rs[0];p<re[0];++p) {
				sub_idx[0] = p;
				for(unsigned int q=rs[1];q<re[1];++q) {
					sub_idx[1] = q;
					for(unsigned int r=rs[2];r<re[2];++r) {
						sub_idx[2] = r;

						upwind.dx(a,b,src_node,sub_idx);
						upwind.dy(c,d,src_node,sub_idx);
						upwind.dz(e,f,src_node,sub_idx);

						CLIP(a,ap,am);	CLIP(b,bp,bm);
						CLIP(c,cp,cm);	CLIP(d,dp,dm);
						CLIP(e,ep,em);	CLIP(f,fp,fm);

						//sign = sign_sussman94(src_node(p,q,r), h2);
						//sign = sign_peng99(src_node(p,q,r), h2);		// Is it possible ?

						const T _phi = (src_node(p,q,r) == 0) ? - eps : src_node(p,q,r);
						const T denom2 = SQ(_phi) + h2 * ( SQ( ((abs(a) > abs(b)) ? a : b) ) 
														 + SQ( ((abs(c) > abs(d)) ? c : d) )
														 + SQ( ((abs(e) > abs(f)) ? e : f) ) );
						sign =  T(_phi / sqrt(denom2));

						if(sign > 0)	{ dst_node(p,q,r) = src_node(p,q,r) - tau * sign * (sqrt(ap*ap + bm*bm + cp*cp + dm*dm + ep*ep + fm*fm) - 1); if(dst_node(p,q,r) < 0) dst_node(p,q,r) = + eps; }
						else			{ dst_node(p,q,r) = src_node(p,q,r) - tau * sign * (sqrt(am*am + bp*bp + cm*cm + dp*dp + em*em + fp*fp) - 1); if(dst_node(p,q,r) > 0) dst_node(p,q,r) = - eps; }
					}					
				}
			}

#undef CLIP
		}

		swap(src, dst);		
		src->synchronize_grid_values();	

	}		// iter

	srl.srl = src; 
	srl.srl0 = dst;
}

};			// fluidspace