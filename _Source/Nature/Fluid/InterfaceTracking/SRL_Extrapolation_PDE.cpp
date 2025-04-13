#include <stdafx.h>
#include <Basis/Basis/Array.h>
#include <Nature/Fluid/InterfaceTracking/SRL_Extrapolation.h>

using namespace std;
using namespace basis;
using namespace fluidspace;

void SRL_Extrapolation_PDE::extrapolate(InterfaceTracking_SRL2D& srl, const Grid2r& _phi, const real zero_value) const {

	const int N = int(srl.srl->n_quadnode);
	const int N_half = N / 2;

	const vector<SRLNode2D<real>>&		nodes = srl.srl->nodes;
	const Grid2<unsigned int, real>&	indices = srl.srl->indices;	

	vector<size_t> internal_indices;

	for(std::size_t i=0;i<nodes.size();++i) {
		index_type ii = nodes[i].index[0], jj = nodes[i].index[1];

		if(_phi(ii,jj) > 0 && _phi(ii,jj+1) > 0 && _phi(ii+1,jj) > 0 && _phi(ii+1,jj+1) > 0) {
			std::cout << ii << " " << jj << " " << std::endl;
			throwError("Error in SRL_Extrapolation_PDE::extrapolate()");
		}

		if(_phi(ii,jj) <= 0 && _phi(ii,jj+1) <= 0 && _phi(ii+1,jj) <= 0 && _phi(ii+1,jj+1) <= 0)
			continue;

		internal_indices.push_back(i);		
	}

	Array<SRLNode2D<real>> phi_nodes(internal_indices.size()), nodes0(internal_indices.size()), nodes1(internal_indices.size());

{

#ifdef USE_OPENMP	
	int my_end = int(internal_indices.size());
	int i;

#pragma omp parallel for
#else	
	std::size_t my_end = internal_indices.size();
	std::size_t i;
#endif
	for( i = 0; i < internal_indices.size(); ++i ) {			
		phi_nodes[i].index[0] = nodes[internal_indices[i]].index[0];
		phi_nodes[i].index[1] = nodes[internal_indices[i]].index[1];

		nodes0[i] = nodes[internal_indices[i]];
		nodes1[i] = nodes[internal_indices[i]];

		for(int p=0;p<N;++p) {
			for(int q=0;q<N;++q) {
				phi_nodes[i](p,q) = _phi(srl.srl->pos(nodes0[i], p, q));
			}
		}
	}
}

// Interpolation
	//for(i = 0; i < my_end; ++i ) {			
	//	for(int p=0;p<N;++p) {
	//		for(int q=0;q<N;++q) {
	//			index_type ii = nodes0[i].index[0], jj = nodes0[i].index[1];

	//			if		( p < N_half && q < N_half && !isbdry(ii  ,jj  )) continue;
	//			else if	( p < N_half && q > N_half && !isbdry(ii  ,jj+1)) continue;
	//			else if	( p > N_half && q < N_half && !isbdry(ii+1,jj  )) continue;
	//			else if	( p > N_half && q > N_half && !isbdry(ii+1,jj+1)) continue;

	//			nodes0[i](p,q) = srl.srl->phi(srl.srl->pos(nodes[internal_indices[i]], p, q));
	//		}
	//	}
	//}	

// Extrapolate
	const real h = srl.srl->h_sub[0];

	const real tau = 0.5f * h;
	const real eps = h, eps2 = eps*eps;

	const size_t nIter = size_t(srl.srl->h / tau + 2);

	for(grid_index_type iter = 0;iter < nIter; ++iter) {

#define CLIP(x,xp,xm) if((x)>0) { xp = (x); xm = 0; } else { xp = 0; xm = (x); }

#ifdef USE_OPENMP	
		int i;
		int my_end = int(nodes0.size());
#pragma omp parallel for
#else	
		size_t i;
		size_t my_end = nodes0.size();
#endif
		for(i = 0; i < my_end; ++i ) {			

			Vector2<real> normal;
			index_type ii = nodes0[i].index[0], jj = nodes0[i].index[1];

			for(int p=0;p<N;++p) {
				for(int q=0;q<N;++q) {

					if		( p < N_half && q < N_half && _phi(ii  ,jj  ) < 0) continue;
					else if	( p < N_half && q > N_half && _phi(ii  ,jj+1) < 0) continue;
					else if	( p > N_half && q < N_half && _phi(ii+1,jj  ) < 0) continue;
					else if	( p > N_half && q > N_half && _phi(ii+1,jj+1) < 0) continue;

					const real phi = phi_nodes[i].val[p][q];						
				
					const real _sign = 1;

					real a(0),b(0),c(0),d(0),phix0(0),phix1(0),phiy0(0),phiy1(0);

					if(p != 0) {
						a	  = srl.srl->hInv_sub[p-1] * (   nodes0[i].val[p][q] -    nodes0[i].val[p-1][q]);
						phix0 = srl.srl->hInv_sub[p-1] * (phi_nodes[i].val[p][q] - phi_nodes[i].val[p-1][q]);
					}
					else if(ii != 0  && indices(ii-1,jj) != SRL<real>::NOT_REFINE) {
						a	  = srl.srl->hInv_sub[N-2] * (   nodes0[i].val[0][q] -                   nodes[indices(ii-1,jj)].val[N-2][q]);
						phix0 = srl.srl->hInv_sub[N-2] * (phi_nodes[i].val[0][q] - _phi(srl.srl->pos(nodes[indices(ii-1,jj)], N-2, q)));
					}

					if(p != N-1) {
						b	  = srl.srl->hInv_sub[p] * (   nodes0[i].val[p+1][q] -    nodes0[i].val[p][q]);
						phix1 = srl.srl->hInv_sub[p] * (phi_nodes[i].val[p+1][q] - phi_nodes[i].val[p][q]);
					}
					else if(indices(ii+1,jj) != SRL<real>::NOT_REFINE) {
						b	  = srl.srl->hInv_sub[0] * (                  nodes[indices(ii+1,jj)].val[1][q] -  nodes0[i].val[N-1][q]);
						phix1 = srl.srl->hInv_sub[0] * (_phi(srl.srl->pos(nodes[indices(ii+1,jj)], 1, q)) - phi_nodes[i].val[N-1][q]);
					}

					if(q != 0) {
						c	  = srl.srl->hInv_sub[q-1] * (   nodes0[i].val[p][q] -    nodes0[i].val[p][q-1]);
						phiy0 = srl.srl->hInv_sub[q-1] * (phi_nodes[i].val[p][q] - phi_nodes[i].val[p][q-1]);
					}
					else if(jj != 0  && indices(ii,jj-1) != SRL<real>::NOT_REFINE) {
						c	  = srl.srl->hInv_sub[N-2] * (   nodes0[i].val[p][0] -                   nodes[indices(ii,jj-1)].val[p][N-2]);
						phiy0 = srl.srl->hInv_sub[N-2] * (phi_nodes[i].val[p][0] - _phi(srl.srl->pos(nodes[indices(ii,jj-1)], p, N-2)));
					}

					if(q != N-1) {
						d	  = srl.srl->hInv_sub[q] * (   nodes0[i].val[p][q+1] -    nodes0[i].val[p][q]);
						phiy1 = srl.srl->hInv_sub[q] * (phi_nodes[i].val[p][q+1] - phi_nodes[i].val[p][q]);
					}
					else if(indices(ii,jj+1) != SRL<real>::NOT_REFINE) {
						d	  = srl.srl->hInv_sub[0] * (                  nodes[indices(ii,jj+1)].val[p][1] -  nodes0[i].val[p][N-1]);					
						phiy1 = srl.srl->hInv_sub[0] * (_phi(srl.srl->pos(nodes[indices(ii,jj+1)], p, 1)) - phi_nodes[i].val[p][N-1]);
					}

					normal[0] = (::abs(phix0) > ::abs(phix1)) ? phix0 : phix1;
					normal[1] = (::abs(phiy0) > ::abs(phiy1)) ? phiy0 : phiy1;
					normal.normalize();
					
					const real snx = _sign * normal[0];
					const real sny = _sign * normal[1];
					real snxp, snxm, snyp, snym;

					CLIP(snx, snxp, snxm);
					CLIP(sny, snyp, snym);

					nodes1[i].val[p][q] = nodes0[i].val[p][q] - tau * (snxp * a + snxm * b + snyp * c + snym * d);	
				}
			}			

		}

		swap(nodes0, nodes1);		

		{
			// Copy
#ifdef USE_OPENMP
			int i;
			int my_end = int(internal_indices.size());
#pragma omp parallel for
#else	
			size_t i;
			size_t my_end = internal_indices.size();
#endif
			for(i = 0; i < my_end; ++i ) {
				srl.srl->nodes[internal_indices[i]] = nodes0[i];
			}
		}

#undef CLIP
	}	
}

void SRL_Extrapolation_PDE::extrapolate(InterfaceTracking_SRL3D& srl, const Grid3r& _phi, const real zero_value) const {
	const int N = int(srl.srl->n_quadnode);
	const int N_half = N / 2;

	const vector<SRLNode3D<real>>&		nodes = srl.srl->nodes;
	const Grid3<unsigned int, real>&	indices = srl.srl->indices;	

	vector<size_t> internal_indices;

	for(std::size_t i=0;i<nodes.size();++i) {
		index_type ii = nodes[i].index[0], jj = nodes[i].index[1], kk = nodes[i].index[2];

		if(_phi(ii,jj,kk  ) > 0 && _phi(ii,jj+1,kk  ) > 0 && _phi(ii+1,jj,kk  ) > 0 && _phi(ii+1,jj+1,kk  ) > 0
		&& _phi(ii,jj,kk+1) > 0 && _phi(ii,jj+1,kk+1) > 0 && _phi(ii+1,jj,kk+1) > 0 && _phi(ii+1,jj+1,kk+1) > 0 ) {

			std::cout << ii << " " << jj << " " << kk << std::endl;
			throwError("Error in SRL_Extrapolation_PDE::extrapolate()");
		}
			

		if(_phi(ii,jj,kk  ) <= 0 && _phi(ii,jj+1,kk  ) <= 0 && _phi(ii+1,jj,kk  ) <= 0 && _phi(ii+1,jj+1,kk  ) <= 0
		&& _phi(ii,jj,kk+1) <= 0 && _phi(ii,jj+1,kk+1) <= 0 && _phi(ii+1,jj,kk+1) <= 0 && _phi(ii+1,jj+1,kk+1) <= 0 )
			continue;

		internal_indices.push_back(i);		
	}

	Array<SRLNode3D<real>> phi_nodes(internal_indices.size()), nodes0(internal_indices.size()), nodes1(internal_indices.size());

	{

#ifdef USE_OPENMP	
		int my_end = int(internal_indices.size());
		int i;

#pragma omp parallel for
#else	
		std::size_t my_end = internal_indices.size();
		std::size_t i;
#endif
		for( i = 0; i < internal_indices.size(); ++i ) {			
			phi_nodes[i].index[0] = nodes[internal_indices[i]].index[0];
			phi_nodes[i].index[1] = nodes[internal_indices[i]].index[1];
			phi_nodes[i].index[2] = nodes[internal_indices[i]].index[2];

			nodes0[i] = nodes[internal_indices[i]];
			nodes1[i] = nodes[internal_indices[i]];

			for(int p=0;p<N;++p) 
				for(int q=0;q<N;++q) 
					for(int r=0;r<N;++r) 
						phi_nodes[i](p,q,r) = _phi(srl.srl->pos(nodes0[i], p, q, r));
		}
	}

	// Extrapolate
	const real h = srl.srl->h_sub[0];

	const real tau = 0.5f * h;
	const real eps = h, eps2 = eps*eps;

	const size_t nIter = size_t(srl.srl->h / tau + 2);

	for(grid_index_type iter = 0;iter < nIter; ++iter) {

#define CLIP(x,xp,xm) if((x)>0) { xp = (x); xm = 0; } else { xp = 0; xm = (x); }

#ifdef USE_OPENMP	
		int i;
		int my_end = int(nodes0.size());
#pragma omp parallel for
#else	
		size_t i;
		size_t my_end = nodes0.size();
#endif
		for(i = 0; i < my_end; ++i ) {			

			Vector3<real> normal;
			index_type ii = nodes0[i].index[0], jj = nodes0[i].index[1], kk = nodes0[i].index[2];

			for(int p=0;p<N;++p) {
				for(int q=0;q<N;++q) {
					for(int r=0;r<N;++r) {

						if		( p < N_half && q < N_half && r < N_half && _phi(ii  ,jj  ,kk  ) < 0) continue;
						else if	( p < N_half && q > N_half && r < N_half && _phi(ii  ,jj+1,kk  ) < 0) continue;
						else if	( p > N_half && q < N_half && r < N_half && _phi(ii+1,jj  ,kk  ) < 0) continue;
						else if	( p > N_half && q > N_half && r < N_half && _phi(ii+1,jj+1,kk  ) < 0) continue;

						else if	( p < N_half && q < N_half && r > N_half && _phi(ii  ,jj  ,kk+1) < 0) continue;
						else if	( p < N_half && q > N_half && r > N_half && _phi(ii  ,jj+1,kk+1) < 0) continue;
						else if	( p > N_half && q < N_half && r > N_half && _phi(ii+1,jj  ,kk+1) < 0) continue;
						else if	( p > N_half && q > N_half && r > N_half && _phi(ii+1,jj+1,kk+1) < 0) continue;

						const real phi = phi_nodes[i].val[p][q][r];						

						const real _sign = 1;

						real a(0),b(0),c(0),d(0),e(0),f(0), phix0(0),phix1(0),phiy0(0),phiy1(0),phiz0(0),phiz1(0);

						if(p != 0) {
							a	  = srl.srl->hInv_sub[p-1] * (   nodes0[i].val[p][q][r] -    nodes0[i].val[p-1][q][r]);
							phix0 = srl.srl->hInv_sub[p-1] * (phi_nodes[i].val[p][q][r] - phi_nodes[i].val[p-1][q][r]);
						}
						else if(ii != 0  && indices(ii-1,jj,kk) != SRL<real>::NOT_REFINE) {
							a	  = srl.srl->hInv_sub[N-2] * (   nodes0[i].val[0][q][r] -                   nodes[indices(ii-1,jj,kk)].val[N-2][q][r]);
							phix0 = srl.srl->hInv_sub[N-2] * (phi_nodes[i].val[0][q][r] - _phi(srl.srl->pos(nodes[indices(ii-1,jj,kk)],    N-2, q, r)));
						}

						if(p != N-1) {
							b	  = srl.srl->hInv_sub[p] * (   nodes0[i].val[p+1][q][r] -    nodes0[i].val[p][q][r]);
							phix1 = srl.srl->hInv_sub[p] * (phi_nodes[i].val[p+1][q][r] - phi_nodes[i].val[p][q][r]);
						}
						else if(indices(ii+1,jj,kk) != SRL<real>::NOT_REFINE) {
							b	  = srl.srl->hInv_sub[0] * (                  nodes[indices(ii+1,jj,kk)].val[1][q][r] -  nodes0[i].val[N-1][q][r]);
							phix1 = srl.srl->hInv_sub[0] * (_phi(srl.srl->pos(nodes[indices(ii+1,jj,kk)], 1, q, r)) - phi_nodes[i].val[N-1][q][r]);
						}

						if(q != 0) {
							c	  = srl.srl->hInv_sub[q-1] * (   nodes0[i].val[p][q][r] -    nodes0[i].val[p][q-1][r]);
							phiy0 = srl.srl->hInv_sub[q-1] * (phi_nodes[i].val[p][q][r] - phi_nodes[i].val[p][q-1][r]);
						}
						else if(jj != 0  && indices(ii,jj-1,kk) != SRL<real>::NOT_REFINE) {
							c	  = srl.srl->hInv_sub[N-2] * (   nodes0[i].val[p][0][r] -                   nodes[indices(ii,jj-1,kk)].val[p][N-2][r]);
							phiy0 = srl.srl->hInv_sub[N-2] * (phi_nodes[i].val[p][0][r] - _phi(srl.srl->pos(nodes[indices(ii,jj-1,kk)],    p, N-2, r)));
						}

						if(q != N-1) {
							d	  = srl.srl->hInv_sub[q] * (   nodes0[i].val[p][q+1][r] -    nodes0[i].val[p][q][r]);
							phiy1 = srl.srl->hInv_sub[q] * (phi_nodes[i].val[p][q+1][r] - phi_nodes[i].val[p][q][r]);
						}
						else if(indices(ii,jj+1,kk) != SRL<real>::NOT_REFINE) {
							d	  = srl.srl->hInv_sub[0] * (                  nodes[indices(ii,jj+1,kk)].val[p][1][r] -  nodes0[i].val[p][N-1][r]);					
							phiy1 = srl.srl->hInv_sub[0] * (_phi(srl.srl->pos(nodes[indices(ii,jj+1,kk)], p, 1, r)) - phi_nodes[i].val[p][N-1][r]);
						}

						if(r != 0) {
							e	  = srl.srl->hInv_sub[r-1] * (   nodes0[i].val[p][q][r] -    nodes0[i].val[p][q][r-1]);
							phiz0 = srl.srl->hInv_sub[r-1] * (phi_nodes[i].val[p][q][r] - phi_nodes[i].val[p][q][r-1]);
						}
						else if(kk != 0  && indices(ii,jj,kk-1) != SRL<real>::NOT_REFINE) {
							e	  = srl.srl->hInv_sub[N-2] * (   nodes0[i].val[p][q][0] -                   nodes[indices(ii,jj,kk-1)].val[p][q][N-2]);
							phiz0 = srl.srl->hInv_sub[N-2] * (phi_nodes[i].val[p][q][0] - _phi(srl.srl->pos(nodes[indices(ii,jj,kk-1)],    p, q, N-2)));
						}

						if(r != N-1) {
							f	  = srl.srl->hInv_sub[r] * (   nodes0[i].val[p][q][r+1] -    nodes0[i].val[p][q][r]);
							phiz1 = srl.srl->hInv_sub[r] * (phi_nodes[i].val[p][q][r+1] - phi_nodes[i].val[p][q][r]);
						}
						else if(indices(ii,jj,kk+1) != SRL<real>::NOT_REFINE) {
							f	  = srl.srl->hInv_sub[0] * (                  nodes[indices(ii,jj,kk+1)].val[p][q][1] -  nodes0[i].val[p][q][N-1]);					
							phiz1 = srl.srl->hInv_sub[0] * (_phi(srl.srl->pos(nodes[indices(ii,jj,kk+1)], p, q, 1)) - phi_nodes[i].val[p][q][N-1]);
						}

						normal[0] = (::abs(phix0) > ::abs(phix1)) ? phix0 : phix1;
						normal[1] = (::abs(phiy0) > ::abs(phiy1)) ? phiy0 : phiy1;
						normal[2] = (::abs(phiz0) > ::abs(phiz1)) ? phiz0 : phiz1;
						normal.normalize();

						const real snx = _sign * normal[0];
						const real sny = _sign * normal[1];
						const real snz = _sign * normal[2];
						real snxp, snxm, snyp, snym, snzp, snzm;

						CLIP(snx, snxp, snxm);
						CLIP(sny, snyp, snym);
						CLIP(snz, snzp, snzm);

						nodes1[i].val[p][q][r] = nodes0[i].val[p][q][r] - tau * (snxp * a + snxm * b + snyp * c + snym * d + snzp * e + snzm * f);	

					}					
				}
			}			

		}

		swap(nodes0, nodes1);		

		{
			// Copy
#ifdef USE_OPENMP
			int i;
			int my_end = int(internal_indices.size());
#pragma omp parallel for
#else	
			size_t i;
			size_t my_end = internal_indices.size();
#endif
			for(i = 0; i < my_end; ++i ) {
				srl.srl->nodes[internal_indices[i]] = nodes0[i];
			}
		}

#undef CLIP
	}
}

