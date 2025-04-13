#include <stdafx.h>
#include <Nature/Fluid/InterfaceTracking/SRL_Redistancing.h>

using namespace std;
using namespace fluidspace;

void SRL_Redistancing_Direct::redistance(InterfaceTracking_SRL2D& srl) const {

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

	// Re-Distancing _ Extract Contour
	vector<Segment2D<real>> tmp;
	vector<vector<Segment2D<real>>> segments;
	segments.resize(srl.srl->nodes.size());

	Vector2<real> org, h;

	real _phi[4];
	bool bdry_xm, bdry_xp, bdry_ym, bdry_yp;

	Grid2<bool,T> isbdry; isbdry.resize(srl.srl->phi);
	isbdry.set(false);
	isbdry.set_boundary(true);

	for(index_type i=0;i<srl.srl->phi.nX-1;++i)	{
		for(index_type j=0;j<srl.srl->phi.nY-1;++j)	{

			if(isbdry(i,j) && isbdry(i+1,j) && isbdry(i,j+1) && isbdry(i+1,j+1))
				continue;

			if(srl.srl->indices(i,j) == SRL2D<real>::NOT_REFINE) {					
				continue;
			}

			const int N = srl.srl->n_quadnode;
			const int half_N = N/2;

			const SRLNode2D<real>& node = srl.srl->nodes[srl.srl->indices(i,j)];

			for(int p=0;p<N-1;++p) {
				for(int q=0;q<N-1;++q) {

					if		( p <  half_N && q <  half_N && isbdry(i  ,j  ) ) continue;
					else if	( p >= half_N && q <  half_N && isbdry(i+1,j  ) ) continue;
					else if	( p <  half_N && q >= half_N && isbdry(i  ,j+1) ) continue;
					else if	( p >= half_N && q >= half_N && isbdry(i+1,j+1) ) continue;

					int i_ = (p < half_N) ? 0 : 1;
					int j_ = (q < half_N) ? 0 : 1;

					bdry_xm = bdry_xp = bdry_ym = bdry_yp = false;

					if(p==half_N-1)		bdry_xp = isbdry(i+1,j+j_);
					else if(p==half_N)	bdry_xm = isbdry(i  ,j+j_);

					if(q==half_N-1)		bdry_yp = isbdry(i+i_,j+1);
					else if(q==half_N)	bdry_ym = isbdry(i+i_,j  );

					org = srl.srl->pos(node, p, q);
					h[0] = srl.srl->h_sub[p]; h[1] = srl.srl->h_sub[q];

					_phi[0] = node(p  ,q  );	_phi[1] = node(p+1,q  );
					_phi[2] = node(p+1,q+1);	_phi[3] = node(p  ,q+1);

					//contouring(tmp, _phi, org, h, bdry_xm, bdry_xp, bdry_ym, bdry_yp);					
					contouring(tmp, _phi, org, h, false, false, false, false);					

					vectorLibrary::append(segments[srl.srl->indices(i,j)], tmp);
				}
			}
		}
	}	

	for(index_type i=0;i<dst->phi.nX;++i)	{
		for(index_type j=0;j<dst->phi.nY;++j)	{
			dst->phi(i,j) = SIGN(dst->phi(i,j)) * (dst->narrowbandWidth);
		}}

	// Re-Distancing _ Measure Distance Directly
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

		const SRLNode2D<T>&	src_node = src->nodes[i];
		SRLNode2D<T>&		dst_node = dst->nodes[i];	

		Vector2<unsigned int> rs, re;
		src->get_updating_subgrid_index(rs, re, src_node.index);

		for(unsigned int p=rs[0];p<re[0];++p) 
			for(unsigned int q=rs[1];q<re[1];++q) 
				dst_node.val[p][q] = SIGN(src_node.val[p][q]) * get_distance(srl.srl->pos(dst_node, p, q), *src, segments);
	}

	dst->synchronize_grid_values();	

	srl.srl = dst; 
	srl.srl0 = src;
}

const real SRL_Redistancing_Direct::get_distance(const Vector2<real>& p, const SRL2D<real>& srl, const vector<vector<Segment2D<real>>>& segments_all) const {

	const index_type _band = index_type(band);
	const index_type nX = srl.phi.nX, nY = srl.phi.nY;

	index_type ii, jj;
	srl.phi.getClosestIdx(ii, jj, p);

	real min_dist = 1.0e+8;

	for(index_type i=MAX(ii-2*_band, 0);i<MIN(ii+2*_band, nX);++i) {
	for(index_type j=MAX(jj-2*_band, 0);j<MIN(jj+2*_band, nY);++j) {

		if(srl.indices(i,j) == SRL2D<real>::NOT_REFINE)
			continue;
	
		const vector<Segment2D<real>>& segments = segments_all[srl.indices(i,j)];

		for(size_t k=0;k<segments.size();++k) {
			const real sqdist = DistanceLibrary::SqDistPointLine(p, segments[k].v[0], segments[k].v[1]);
			min_dist = MIN(min_dist, sqrt(sqdist));
		}

	}}

	return min_dist;
}	

void SRL_Redistancing_Direct::redistance(InterfaceTracking_SRL3D& srl) const {

	typedef real T;
	typedef TriMesh_Simple<Vector3<T>, vector<Vector3<T>>> Mesh;

	SRL3D<real>* src = srl.srl, *dst = srl.srl0;

	// Re-distancing Coarse Grid
	src->modify_coarse_phi_except_near_interface();
	geometryspace::Levelset_Redistancing_PDE<real>::redistance_pde(dst->phi, src->phi, Redistancing_EE<Upwind_1st<Grid3<T,T>>>(*this));

	// Extract Narrow-band Cells using Coarse Grid Phi
	SRL3D<T>::extract_narrowband(dst->indices, dst->phi, src->narrowbandWidth, unsigned int(SRL<T>::NEED_REFINE), unsigned int(SRL<T>::NOT_REFINE));

	// Add Narrow-band Cells which have sub-cell interfaces
	srl.add_narrowband_subcell_interface(dst->indices, src->nodes, unsigned int(SRL<T>::NEED_REFINE));

	// Subtract Narrow-band Cells which is contained in Kinematic Boundary
	srl.subtract_narrowband_kinematicbdry(srl.srl0->indices, unsigned int(SRL<real>::NOT_REFINE));

	// Refine Grid
	dst->refine_subgrid(dst->indices, unsigned int(SRL<T>::NEED_REFINE), unsigned int(SRL<T>::NOT_REFINE));

	// Copy Values & Structures
	dst->set_value(*src);	dst->synchronize_grid_values();
	src->set(*dst);

	// Marching Interface
	src->createRenderingPrimitive(Grid3<bool,real>());

	AABB<T> bbox = src->mesh.getAABB();
	bbox.expand(bbox.extent() * T(0.1));
	
	Octree_Hashing<int, Octree_Hashing_AddProcess_TriMesh<Mesh>, T>* octree = new Octree_Hashing<int, Octree_Hashing_AddProcess_TriMesh<Mesh>, T>(bbox, 9);

	Octree_Hashing_AddProcess_TriMesh<Mesh> addProcess(src->mesh);
	for(unsigned int i=0;i<src->mesh.numFaces;++i)
		octree->add(int(i), addProcess);

	for(index_type i=0;i<dst->phi.nX;++i)	{
		for(index_type j=0;j<dst->phi.nY;++j)	{
			for(index_type k=0;k<dst->phi.nZ;++k)	{
				dst->phi(i,j,k) = SIGN(dst->phi(i,j,k)) * (dst->narrowbandWidth);
			}}}

	// Re-Distancing _ Measure Distance Directly
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

		const SRLNode3D<T>&	src_node = src->nodes[i];
		SRLNode3D<T>&		dst_node = dst->nodes[i];	

		Vector3<unsigned int> rs, re;
		src->get_updating_subgrid_index(rs, re, src_node.index);

		for(unsigned int p=rs[0];p<re[0];++p) {
			for(unsigned int q=rs[1];q<re[1];++q) {
				for(unsigned int r=rs[2];r<re[2];++r) {

					const Vector3<T> pos = srl.srl->pos(dst_node, p, q, r);

					Octree_Hashing_MinDistProcess_TriMesh<Mesh> min_dist(src->mesh);
					octree->lookup(pos, min_dist);

					dst_node.val[p][q][r] = SIGN(src_node.val[p][q][r]) * sqrt(min_dist.minSqDist);
				}}}
	}

	dst->synchronize_grid_values();	

	srl.srl = dst; 
	srl.srl0 = src;

	src->mesh.release();

	delete octree;
}
