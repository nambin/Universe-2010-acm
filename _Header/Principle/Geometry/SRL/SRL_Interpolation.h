#pragma once

namespace geometryspace {

template<typename T, int N>
const T SRL<T,N>::SRL_Interpolator::operator()(const SRL2D<T,N>& srl, index_type i, index_type j, const Vector2<T>& x) const {
	num_itrp++;

	T w[2][10], vy, ret;

	srl.get_weight(w[0], x[0]);
	srl.get_weight(w[1], x[1]);

	const SRLNode2D<T,N>& node = srl.node(i,j);

	ret = 0;
	for(int p=0;p<srl.n_quadnode;++p) {
		vy = 0;
		for(int q=0;q<srl.n_quadnode;++q) {
			vy += (node.val[p][q] * w[1][q]);
		}
		ret += (vy * w[0][p]);
	}

	if(clipping) {
		const int p = srl.get_sub_index(x[0]);
		const int q = srl.get_sub_index(x[1]);

		const SRLNode2D<T,N>& node = srl.node(i,j);

		if((ret > node(p,q) && ret > node(p,q+1) && ret > node(p+1,q) && ret > node(p+1,q+1))
		|| (ret < node(p,q) && ret < node(p,q+1) && ret < node(p+1,q) && ret < node(p+1,q+1))) {

			num_clipping++;

			T w[2];
			w[0] = (x[0] - srl.r_sub[p]) * srl.hInv_norm_sub[p];
			w[1] = (x[1] - srl.r_sub[q]) * srl.hInv_norm_sub[q];

			return (1-w[0]) * ((1-w[1])*node(p  ,q)	+ w[1]*node(p  ,q+1))
					+ w[0]  * ((1-w[1])*node(p+1,q)	+ w[1]*node(p+1,q+1)); 
		}
		else
			return ret;
	}

	return ret;
}

template<typename T, int N>
const T SRL<T,N>::SRL_Lerp::operator()(const SRL2D<T,N>& srl, index_type i, index_type j, const Vector2<T>& x) const {
	T w[2];

	const int p = srl.get_sub_index(x[0]);
	const int q = srl.get_sub_index(x[1]);

	const SRLNode2D<T,N>& node = srl.node(i,j);

	w[0] = (x[0] - srl.r_sub[p]) * srl.hInv_norm_sub[p];
	w[1] = (x[1] - srl.r_sub[q]) * srl.hInv_norm_sub[q];

	return (1-w[0]) * ((1-w[1])*node(p  ,q)	+ w[1]*node(p  ,q+1))
			+ w[0]  * ((1-w[1])*node(p+1,q)	+ w[1]*node(p+1,q+1)); 
}

template<typename T, int N>
const T SRL<T,N>::SRL_Coarse_Lerp::operator()(const SRL2D<T,N>& srl, index_type i, index_type j, const Vector2<T>& x) const {
	return srl.phi.bilerp(i,j,x[0],x[1]);
}

template<typename T, int N>
const T SRL<T,N>::SRL_Interpolator::operator()(const SRL3D<T,N>& srl, index_type i, index_type j, index_type k, const Vector3<T>& x) const {
	num_itrp++;

	T w[3][10], vy, vz, ret;

	srl.get_weight(w[0], x[0]);
	srl.get_weight(w[1], x[1]);
	srl.get_weight(w[2], x[2]);

	const SRLNode3D<T,N>& node = srl.node(i,j,k);

	ret = 0;
	for(int p=0;p<srl.n_quadnode;++p) {
		vy = 0;
		for(int q=0;q<srl.n_quadnode;++q) {
			vz = 0;
			for(int r=0;r<srl.n_quadnode;++r) {
				vz += (node.val[p][q][r] * w[2][r]);
			}
			vy += (vz * w[1][q]);
		}
		ret += (vy * w[0][p]);
	}

	if(clipping) {
		const int p = srl.get_sub_index(x[0]);
		const int q = srl.get_sub_index(x[1]);
		const int r = srl.get_sub_index(x[2]);

		const SRLNode3D<T,N>& node = srl.node(i,j,k);

		if((ret > node(p,q,r  ) && ret > node(p,q+1,r  ) && ret > node(p+1,q,r  ) && ret > node(p+1,q+1,r  )
		 && ret > node(p,q,r+1) && ret > node(p,q+1,r+1) && ret > node(p+1,q,r+1) && ret > node(p+1,q+1,r+1) )
		 || 
		   (ret < node(p,q,r  ) && ret < node(p,q+1,r  ) && ret < node(p+1,q,r  ) && ret < node(p+1,q+1,r  )
	     && ret < node(p,q,r+1) && ret < node(p,q+1,r+1) && ret < node(p+1,q,r+1) && ret < node(p+1,q+1,r+1) )
		  ) 
		{
				num_clipping++;

				T w[3];
				w[0] = (x[0] - srl.r_sub[p]) * srl.hInv_norm_sub[p];
				w[1] = (x[1] - srl.r_sub[q]) * srl.hInv_norm_sub[q];
				w[2] = (x[2] - srl.r_sub[r]) * srl.hInv_norm_sub[r];

				return (1-w[0]) * ((1-w[1])* ((1-w[2])*node(p  ,q  ,r) + w[2]*node(p  ,q  ,r+1))	
									+ w[1] * ((1-w[2])*node(p  ,q+1,r) + w[2]*node(p  ,q+1,r+1)) ) 
					    + w[0]  * ((1-w[1])* ((1-w[2])*node(p+1,q  ,r) + w[2]*node(p+1,q  ,r+1))
									+ w[1] * ((1-w[2])*node(p+1,q+1,r) + w[2]*node(p+1,q+1,r+1)) ); 
		}
		else
			return ret;
	}

	return ret;
}

template<typename T, int N>
const T SRL<T,N>::SRL_Lerp::operator()(const SRL3D<T,N>& srl, index_type i, index_type j, index_type k, const Vector3<T>& x) const {
	T w[3];

	const int p = srl.get_sub_index(x[0]);
	const int q = srl.get_sub_index(x[1]);
	const int r = srl.get_sub_index(x[2]);

	const SRLNode3D<T,N>& node = srl.node(i,j,k);

	w[0] = (x[0] - srl.r_sub[p]) * srl.hInv_norm_sub[p];
	w[1] = (x[1] - srl.r_sub[q]) * srl.hInv_norm_sub[q];
	w[2] = (x[2] - srl.r_sub[r]) * srl.hInv_norm_sub[r];

	return (1-w[0]) * ((1-w[1])* ((1-w[2])*node(p  ,q  ,r) + w[2]*node(p  ,q  ,r+1))	
						+ w[1] * ((1-w[2])*node(p  ,q+1,r) + w[2]*node(p  ,q+1,r+1)) ) 
			+ w[0]  * ((1-w[1])* ((1-w[2])*node(p+1,q  ,r) + w[2]*node(p+1,q  ,r+1))
						+ w[1] * ((1-w[2])*node(p+1,q+1,r) + w[2]*node(p+1,q+1,r+1)) ); 
}

template<typename T, int N> template<typename T2>
const T2 SRL<T,N>::SRL_Lerp::operator()(const SRL3D<T,N>& srl, const SRLNode3D<T2,N>& node, const Vector3<T>& x) const {
	T w[3];

	const int p = srl.get_sub_index(x[0]);
	const int q = srl.get_sub_index(x[1]);
	const int r = srl.get_sub_index(x[2]);

	w[0] = (x[0] - srl.r_sub[p]) * srl.hInv_norm_sub[p];
	w[1] = (x[1] - srl.r_sub[q]) * srl.hInv_norm_sub[q];
	w[2] = (x[2] - srl.r_sub[r]) * srl.hInv_norm_sub[r];

	return (1-w[0]) * ((1-w[1])* ((1-w[2])*node(p  ,q  ,r) + w[2]*node(p  ,q  ,r+1))	
						+ w[1] * ((1-w[2])*node(p  ,q+1,r) + w[2]*node(p  ,q+1,r+1)) ) 
			+ w[0]  * ((1-w[1])* ((1-w[2])*node(p+1,q  ,r) + w[2]*node(p+1,q  ,r+1))
						+ w[1] * ((1-w[2])*node(p+1,q+1,r) + w[2]*node(p+1,q+1,r+1)) ); 
}

template<typename T, int N>
const T SRL<T,N>::SRL_Coarse_Lerp::operator()(const SRL3D<T,N>& srl, index_type i, index_type j, index_type k, const Vector3<T>& x) const {
	return srl.phi.trilerp(i,j,k,x[0],x[1],x[2]);
}

template<typename T, int N> template<typename T2>
const T2 SRL<T,N>::SRL_Coarse_Lerp::operator()(const Grid3<T2,T>& grid, index_type i, index_type j, index_type k, const Vector3<T>& x) const {
	return grid.trilerp(i,j,k,x[0],x[1],x[2]);
}

};		// geometryspace