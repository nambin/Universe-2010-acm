#pragma once

#include <hash_map>
#include "SRL_FileIO.h"
#include "SRL_Derivative.h"

namespace geometryspace {

template<typename T, int N>
void SRL3D<T,N>::getNormals_SubGrid(std::vector<SRLNode3D<Vector3<T>>>& normals) {
	
	normals.resize(nodes.size());

#ifdef USE_OPENMP	
	int my_end = int(nodes.size());
	int k;

#pragma omp parallel for
#else	
	size_t my_end = nodes.size();
	size_t k;
#endif
	for(k=0;k<my_end;++k) {	

		Derivative_SRL_2nd<SRL3D<real>> dev(*this);

		index_type i[3];
		unsigned int index[3];
		Vector3<unsigned int> rs, re;

		for(size_t x=0;x<3;++x) {
			i[x] = nodes[k].index[x];
			normals[k].index[x] = nodes[k].index[x];
		}

		rs[0] = rs[1] = rs[2] = 0;
		re[0] = re[1] = re[2] = N;

		if(i[0] == 0)				rs[0] = N / 2;
		else if(i[0] == phi.nX-2)	re[0] = N / 2 + 1;

		if(i[1] == 0)				rs[1] = N / 2;
		else if(i[1] == phi.nY-2)	re[1] = N / 2 + 1;

		if(i[2] == 0)				rs[2] = N / 2;
		else if(i[2] == phi.nZ-2)	re[2] = N / 2 + 1;

		for(unsigned int p=rs[0];p<re[0];++p) {
			index[0] = p;
			for(unsigned int q=rs[1];q<re[1];++q) {
				index[1] = q;
				for(unsigned int r=rs[2];r<re[2];++r) {
					index[2] = r;

					normals[k].val[p][q][r] = unitVector(Vector3<T>(dev.dx(nodes[k], index), dev.dy(nodes[k], index), dev.dz(nodes[k], index)));
				}
			}
		}

		if( i[0] == 0 )				normals[k].neumann_x_lower();
		else if( i[0] == phi.nX-2 )	normals[k].neumann_x_upper();

		if( i[1] == 0 )				normals[k].neumann_y_lower();
		else if( i[1] == phi.nY-2 )	normals[k].neumann_y_upper();

		if( i[2] == 0 )				normals[k].neumann_z_lower();
		else if( i[2] == phi.nZ-2 )	normals[k].neumann_z_upper();
	}
}

template<typename T, int N>
const T SRL3D<T,N>::getTotalVolume() const {
	const T zero(0.0), one(1.0);

	T vol = 0;
	for(index_type i=1;i<phi.nX-1;++i)
		for(index_type j=1;j<phi.nY-1;++j)
			for(index_type k=1;k<phi.nZ-1;++k)
				vol += MacGridLibrary<T>::heaviside(phi(i,j,k) * hInv, zero, one);;

	return vol * h * h * h;
}

template<typename T, int N>
const T SRL3D<T,N>::getTotalVolume(const Grid3<bool,T>& isbdry) const {
	const T zero(0.0), one(1.0);

	T vol = 0;
	for(index_type i=1;i<phi.nX-1;++i)
		for(index_type j=1;j<phi.nY-1;++j)
			for(index_type k=1;k<phi.nZ-1;++k)
				if(isbdry(i,j,k) == false)
					vol += MacGridLibrary<T>::heaviside(phi(i,j,k) * hInv, zero, one);;

	return vol * h * h * h;
}

template<typename T, int N>
void SRL3D<T,N>::createRenderingPrimitive() {
	SRL<real>::SRL_Interpolator interp(true);
	SRL<real>::SRL_Lerp			lerp;
	SRL<real>::SRL_Coarse_Lerp	coarse_lerp;

	if(nX * nY * nZ < 66 * 66 * 66) {
//	if(nX * nY * nZ < 130 * 130 * 130) {
//	if(nX * nY * nZ < 2) {
		Grid3<T,T> phiNode(nX * (N-1) + 1, nY * (N-1) + 1, nZ * (N-1) + 1,  h / T(N-1), h / T(N-1), h / T(N-1), Vector3<T>(0));
		for(index_type i=0;i<phiNode.nX;++i) {
			for(index_type j=0;j<phiNode.nY;++j) {
				for(index_type k=0;k<phiNode.nZ;++k) {
					index_type ii, jj, kk;	Vector3<T> x;
					get_idx_normalized_coord(ii, jj, kk, x, phiNode.pos(i,j,k));

					if(is_exist_node(ii,jj,kk))	{ phiNode(i,j,k) = lerp(*this, ii, jj, kk, x);			}
					else						{ phiNode(i,j,k) = coarse_lerp(*this, ii, jj, kk, x);	}

				}			
			}
		}

		Levelset3D<T>::marchingCube(mesh, phiNode, Grid3<bool,T>());
	}
	else {
		Grid3<T,T> phiNode(nX+1, nY+1, nZ+1, h, h, h, Vector3<T>(0));
		phi.setNeumannBoundary();
		for(index_type i=0;i<phiNode.nX;++i) {
			for(index_type j=0;j<phiNode.nY;++j) {
				for(index_type k=0;k<phiNode.nZ;++k) {
					phiNode(i,j,k) = 0.125f * (phi(i,j,k  ) + phi(i+1,j,k  ) + phi(i,j+1,k  ) + phi(i+1,j+1,k  )
											+  phi(i,j,k+1) + phi(i+1,j,k+1) + phi(i,j+1,k+1) + phi(i+1,j+1,k+1));
				}			
			}
		}

		Levelset3D<T>::marchingCube(mesh, phiNode, Grid3<bool,T>());
	}	
}

template<typename T, int N>
void SRL3D<T,N>::createRenderingPrimitive(const Grid3<bool,T>& isbdry) {
	//createRenderingPrimitive(); return;

	if(!isbdry.is_empty()) {		
		marchingInterface(mesh, isbdry);
	}
	else {
		Grid3<bool,T> _isbdry; _isbdry.resize(phi);
		_isbdry.set(false);
		_isbdry.set_boundary(true);
		
		marchingInterface(mesh, _isbdry);
	}	
}

template<typename T>
class Hashing_SRL_Marching_Wall {
public : 	
	Hashing_SRL_Marching_Wall(const Grid3<T,T>& _phi) : phi(_phi), numVertices(0) { init_state(); }

	void init_state() { for(int i=0;i<4;++i) data[i] = -1; }

	const int get_key(int i, int j, int k, int p, int q, int r) { 
		if(p == 2) { i++; p = 0; }
		if(q == 2) { j++; q = 0; }
		if(r == 2) { k++; r = 0; }

		return ((phi.idx(i,j,k) * 2 + p) * 2 + q) * 2 + r;
	}

	const int operator()(int i) { 
		if(!map.count(data[i]))	return -1;
		else					return map[data[i]];
	}

	const int include(const Vector3<T>& p, typename Triangle3D<T>::TYPE _type) {
		posVerts.push_back(p);
		vert_types.push_back(_type);
		numVertices++;

		return numVertices-1;
	}

	const int include(int i, const Vector3<T>& p, typename Triangle3D<T>::TYPE _type) {
		posVerts.push_back(p);
		vert_types.push_back(_type);

		if(map.count(data[i]) != 0)
			throwError("Error in Hash _ Marching Interface");

		map[data[i]] = numVertices;
		numVertices++;

		return numVertices-1;
	}	

public : 
	const Grid3<T,T>&								phi;

	int												data[4];
	int												numVertices;

	std::vector<Vector3<T>>							posVerts;
	std::vector<typename Triangle3D<T>::TYPE>		vert_types;

	stdext::hash_map<int, int>						map;
};

template<typename T, int N>
class Hashing_SRL_Marching {
public : 	
	Hashing_SRL_Marching(const Grid3<T,T>& _phi) : phi(_phi), numVertices(0) { init_state(); }

	void init_state() { for(int i=0;i<8;++i) data[i] = -1; }

	void set_data(int ith, int i, int j, int k, int p, int q, int r) { 
		if(p == N-1) { i++; p = 0; }
		if(q == N-1) { j++; q = 0; }
		if(r == N-1) { k++; r = 0; }

		data[ith] = ((phi.idx(i,j,k) * (N-1) + p) * (N-1) + q) * (N-1) + r;
	}

	const int operator()(int e0, int e1) { 
		unsigned long long ei = get_key(e0, e1);

		if(!map.count(ei))	return -1;
		else				return map[ei];
	}

	const int include(const Vector3<T>& p, typename Triangle3D<T>::TYPE _type) {
		posVerts.push_back(p);
		vert_types.push_back(_type);
		numVertices++;

		return numVertices-1;
	}

	const int include(int e0, int e1, const Vector3<T>& p, typename Triangle3D<T>::TYPE _type) {
		unsigned long long ei = get_key(e0, e1);

		posVerts.push_back(p);
		vert_types.push_back(_type);

		if(map.count(ei) != 0)
			throwError("Error in Hash _ Marching Interface");

		map[ei] = numVertices;
		numVertices++;

		return numVertices-1;
	}	

private : 
	unsigned long long get_key(int e0, int e1) const {
		int a = data[e0], b = data[e1];
		if(a > b) std::swap(a, b);

		unsigned long long ei = a;
		ei = ei << 32;	ei += b;

		return ei;
	}

public : 
	const Grid3<T,T>&								phi;

	int												data[8];

	int												numVertices;

	std::vector<Vector3<T>>							posVerts;
	std::vector<typename Triangle3D<T>::TYPE>		vert_types;

	stdext::hash_map<unsigned long long, int>		map;
};

template<typename T, int N>
void SRL3D<T,N>::marchingInterface(TriMesh_Simple<Vector3<T>, vector<Vector3<T>>>& mesh, const Grid3<bool,T>& isbdry) {

	Hashing_SRL_Marching<T,N>		hashing(phi);
	Hashing_SRL_Marching_Wall<T>	x_hashing(phi), y_hashing(phi), z_hashing(phi);

	mesh.init();	

	vector<Triangle3D<T>::TYPE>		types;
	vector<Triangle3D<T>>			tmp, x_wall, y_wall, z_wall;
	Vector3<T>						org, h;

	T _phi[8];
	bool bdry_xm, bdry_xp, bdry_ym, bdry_yp, bdry_zm, bdry_zp;

	for(index_type i=0;i<phi.nX-1;++i)	{
		for(index_type j=0;j<phi.nY-1;++j)	{
			for(index_type k=0;k<phi.nZ-1;++k)	{

				if(isbdry(i,j,k  ) && isbdry(i+1,j,k  ) && isbdry(i,j+1,k  ) && isbdry(i+1,j+1,k  )
				&& isbdry(i,j,k+1) && isbdry(i+1,j,k+1) && isbdry(i,j+1,k+1) && isbdry(i+1,j+1,k+1) )
					continue;

				if(indices(i,j,k) == NOT_REFINE || nodes[indices(i,j,k)].has_interface() == false) {

					if( indices(i,j,k) == NOT_REFINE 
					&& !(phi(i,j,k  ) >  0 && phi(i+1,j,k  ) >  0 && phi(i,j+1,k  ) >  0 && phi(i+1,j+1,k  ) >  0 &&
						 phi(i,j,k+1) >  0 && phi(i+1,j,k+1) >  0 && phi(i,j+1,k+1) >  0 && phi(i+1,j+1,k+1) >  0 ) 
					&& !(phi(i,j,k  ) <= 0 && phi(i+1,j,k  ) <= 0 && phi(i,j+1,k  ) <= 0 && phi(i+1,j+1,k  ) <= 0 &&
						 phi(i,j,k+1) <= 0 && phi(i+1,j,k+1) <= 0 && phi(i,j+1,k+1) <= 0 && phi(i+1,j+1,k+1) <= 0 ) ) {

							cout << i << " " << j << " " << k << endl;
							cout << isbdry(i,j,k  ) << " " << isbdry(i+1,j,k  ) << " " << isbdry(i,j+1,k  ) << " " << isbdry(i+1,j+1,k  ) << endl;
							cout << isbdry(i,j,k+1) << " " << isbdry(i+1,j,k+1) << " " << isbdry(i,j+1,k+1) << " " << isbdry(i+1,j+1,k+1) << endl;
							cout << phi(i,j,k  ) << " " << phi(i+1,j,k  ) << " " << phi(i,j+1,k  ) << " " << phi(i+1,j+1,k  ) << endl;
							cout << phi(i,j,k+1) << " " << phi(i+1,j,k+1) << " " << phi(i,j+1,k+1) << " " << phi(i+1,j+1,k+1) << endl;
							throwError("Error in SRL3D<T,N>::marchingInterface()");

					}					

					if(phi(i,j,k) > 0 || (!isbdry(i,j,k  ) && !isbdry(i+1,j,k  ) && !isbdry(i,j+1,k  ) && !isbdry(i+1,j+1,k  )
									   && !isbdry(i,j,k+1) && !isbdry(i+1,j,k+1) && !isbdry(i,j+1,k+1) && !isbdry(i+1,j+1,k+1) )
						)
						continue;

					_phi[0] = _phi[1] = _phi[2] = _phi[3] = _phi[4] = _phi[5] = _phi[6] = _phi[7] = -1;
					h[0] = phi.dx*0.5f; h[1] = phi.dy*0.5f; h[2] = phi.dz*0.5f;

					hashing.init_state();
					for(index_type p=0;p<2;++p) {
						for(index_type q=0;q<2;++q) {
							for(index_type r=0;r<2;++r) {
								if(isbdry(i+p,j+q,k+r))
									continue;

								bdry_xm = bdry_xp = bdry_ym = bdry_yp = bdry_zm = bdry_zp = false;

								if(p==0)	bdry_xp = isbdry(i+1,j+q,k+r);
								else 		bdry_xm = isbdry(i  ,j+q,k+r);

								if(q==0)	bdry_yp = isbdry(i+p,j+1,k+r);
								else		bdry_ym = isbdry(i+p,j  ,k+r);

								if(r==0)	bdry_zp = isbdry(i+p,j+q,k+1);
								else		bdry_zm = isbdry(i+p,j+q,k  );

								org = phi.pos(i,j,k) + Vector3<T>(h[0]*p, h[1]*q, h[2]*r);

								const Vector3<T> x_axis(h[0],0,0);
								const Vector3<T> y_axis(0,h[1],0);
								const Vector3<T> z_axis(0,0,h[2]);

								int keys[8];
								keys[0] = x_hashing.get_key(i,j,k,p,q,r);   keys[1] = x_hashing.get_key(i,j,k,p+1,q,r);	  keys[2] = x_hashing.get_key(i,j,k,p+1,q,r+1);   keys[3] = x_hashing.get_key(i,j,k,p,q,r+1);
								keys[4] = x_hashing.get_key(i,j,k,p,q+1,r); keys[5] = x_hashing.get_key(i,j,k,p+1,q+1,r); keys[6] = x_hashing.get_key(i,j,k,p+1,q+1,r+1); keys[7] = x_hashing.get_key(i,j,k,p,q+1,r+1);

								tmp.clear();

								if(bdry_xm) { 
									x_hashing.data[0] = keys[0]; x_hashing.data[1] = keys[3]; x_hashing.data[2] = keys[7]; x_hashing.data[3] = keys[4];
									
									if(i==0) contouring_wall(x_wall, x_hashing, _phi[0], _phi[3], _phi[7], _phi[4], org  , +z_axis, +y_axis, Triangle3D<T>::BDRY_XM); 
									else	 contouring_bdry(tmp, hashing	  , _phi[0], _phi[3], _phi[7], _phi[4], org  , +z_axis, +y_axis, Triangle3D<T>::BDRY_XM); 
								}
								if(bdry_xp) { 
									x_hashing.data[0] = keys[6]; x_hashing.data[1] = keys[2]; x_hashing.data[2] = keys[1]; x_hashing.data[3] = keys[5];

									if(i==phi.nX-2)	contouring_wall(x_wall, x_hashing, _phi[6], _phi[2], _phi[1], _phi[5], org+h, -y_axis, -z_axis, Triangle3D<T>::BDRY_XP); 
									else			contouring_bdry(tmp, hashing	 , _phi[6], _phi[2], _phi[1], _phi[5], org+h, -y_axis, -z_axis, Triangle3D<T>::BDRY_XP); 
								}

								if(bdry_ym) { 
									y_hashing.data[0] = keys[0]; y_hashing.data[1] = keys[1]; y_hashing.data[2] = keys[2]; y_hashing.data[3] = keys[3];

									if(j==0) contouring_wall(y_wall, y_hashing, _phi[0], _phi[1], _phi[2], _phi[3], org  , +x_axis, +z_axis, Triangle3D<T>::BDRY_YM); 
									else	 contouring_bdry(tmp, hashing	  , _phi[0], _phi[1], _phi[2], _phi[3], org  , +x_axis, +z_axis, Triangle3D<T>::BDRY_YM); 
								}
								if(bdry_yp) { 
									y_hashing.data[0] = keys[6]; y_hashing.data[1] = keys[5]; y_hashing.data[2] = keys[4]; y_hashing.data[3] = keys[7];

									if(j==phi.nY-2)	contouring_wall(y_wall, y_hashing, _phi[6], _phi[5], _phi[4], _phi[7], org+h, -z_axis, -x_axis, Triangle3D<T>::BDRY_YP); 
									else			contouring_bdry(tmp, hashing	 , _phi[6], _phi[5], _phi[4], _phi[7], org+h, -z_axis, -x_axis, Triangle3D<T>::BDRY_YP); 
								}

								if(bdry_zm) { 
									z_hashing.data[0] = keys[0]; z_hashing.data[1] = keys[4]; z_hashing.data[2] = keys[5]; z_hashing.data[3] = keys[1];

									if(k==0) contouring_wall(z_wall, z_hashing, _phi[0], _phi[4], _phi[5], _phi[1], org  , +y_axis, +x_axis, Triangle3D<T>::BDRY_ZM); 
									else	 contouring_bdry(tmp, hashing	  , _phi[0], _phi[4], _phi[5], _phi[1], org  , +y_axis, +x_axis, Triangle3D<T>::BDRY_ZM); 
								}
								if(bdry_zp) { 
									z_hashing.data[0] = keys[6]; z_hashing.data[1] = keys[7]; z_hashing.data[2] = keys[3]; z_hashing.data[3] = keys[2];

									if(k==phi.nZ-2)	 contouring_wall(z_wall, z_hashing, _phi[6], _phi[7], _phi[3], _phi[2], org+h, -x_axis, -y_axis, Triangle3D<T>::BDRY_ZP); 
									else			 contouring_bdry(tmp, hashing	 , _phi[6], _phi[7], _phi[3], _phi[2], org+h, -x_axis, -y_axis, Triangle3D<T>::BDRY_ZP); 
								}	

								//if(bdry_xm) contouring_bdry(tmp, hashing	 , _phi[0], _phi[3], _phi[7], _phi[4], org  , +z_axis, +y_axis, Triangle3D<T>::BDRY_XM); 
								//if(bdry_xp) contouring_bdry(tmp, hashing	 , _phi[6], _phi[2], _phi[1], _phi[5], org+h, -y_axis, -z_axis, Triangle3D<T>::BDRY_XP); 
								//if(bdry_ym) contouring_bdry(tmp, hashing	 , _phi[0], _phi[1], _phi[2], _phi[3], org  , +x_axis, +z_axis, Triangle3D<T>::BDRY_YM); 
								//if(bdry_yp) contouring_bdry(tmp, hashing	 , _phi[6], _phi[5], _phi[4], _phi[7], org+h, -z_axis, -x_axis, Triangle3D<T>::BDRY_YP); 
								//if(bdry_zm) contouring_bdry(tmp, hashing	 , _phi[0], _phi[4], _phi[5], _phi[1], org  , +y_axis, +x_axis, Triangle3D<T>::BDRY_ZM); 
								//if(bdry_zp) contouring_bdry(tmp, hashing	 , _phi[6], _phi[7], _phi[3], _phi[2], org+h, -x_axis, -y_axis, Triangle3D<T>::BDRY_ZP); 

								for(size_t m=0;m<tmp.size();++m) {

									mesh.v2f.push_back(tmp[m].index[0]);
									mesh.v2f.push_back(tmp[m].index[1]);
									mesh.v2f.push_back(tmp[m].index[2]);
								}		
							}							
						}
					}				
					continue;
				}

				const int half_N = N/2;

				const SRLNode3D<T,N>& node = nodes[indices(i,j,k)];

				if(node.has_interface() == false)
					continue;

				for(int p=0;p<N-1;++p) {
					for(int q=0;q<N-1;++q) {
						for(int r=0;r<N-1;++r) {

							if		( p <  half_N && q <  half_N && r <  half_N && isbdry(i  ,j  ,k  ) ) continue;
							else if	( p >= half_N && q <  half_N && r <  half_N && isbdry(i+1,j  ,k  ) ) continue;
							else if	( p <  half_N && q >= half_N && r <  half_N && isbdry(i  ,j+1,k  ) ) continue;
							else if	( p >= half_N && q >= half_N && r <  half_N && isbdry(i+1,j+1,k  ) ) continue;

							else if	( p <  half_N && q <  half_N && r >= half_N && isbdry(i  ,j  ,k+1) ) continue;
							else if	( p >= half_N && q <  half_N && r >= half_N && isbdry(i+1,j  ,k+1) ) continue;
							else if	( p <  half_N && q >= half_N && r >= half_N && isbdry(i  ,j+1,k+1) ) continue;
							else if	( p >= half_N && q >= half_N && r >= half_N && isbdry(i+1,j+1,k+1) ) continue;

							int i_ = (p < half_N) ? 0 : 1;
							int j_ = (q < half_N) ? 0 : 1;
							int k_ = (r < half_N) ? 0 : 1;

							bdry_xm = bdry_xp = bdry_ym = bdry_yp = bdry_zm = bdry_zp = false;

							if(p==half_N-1)		bdry_xp = isbdry(i+1,j+j_,k+k_);
							else if(p==half_N)	bdry_xm = isbdry(i  ,j+j_,k+k_);

							if(q==half_N-1)		bdry_yp = isbdry(i+i_,j+1,k+k_);
							else if(q==half_N)	bdry_ym = isbdry(i+i_,j  ,k+k_);

							if(r==half_N-1)		bdry_zp = isbdry(i+i_,j+j_,k+1);
							else if(r==half_N)	bdry_zm = isbdry(i+i_,j+j_,k  );

							org = pos(node, p, q, r);
							h[0] = h_sub[p]; h[1] = h_sub[q]; h[2] = h_sub[r];

							_phi[0] = node(p,q  ,r);	_phi[1] = node(p+1,q  ,r);	_phi[2] = node(p+1,q  ,r+1);	_phi[3] = node(p,q  ,r+1);
							_phi[4] = node(p,q+1,r);	_phi[5] = node(p+1,q+1,r);	_phi[6] = node(p+1,q+1,r+1);	_phi[7] = node(p,q+1,r+1);	

							hashing.init_state();
							hashing.set_data(0,i,j,k,p,q  ,r); hashing.set_data(1,i,j,k,p+1,q  ,r); hashing.set_data(2,i,j,k,p+1,q  ,r+1); hashing.set_data(3,i,j,k,p,q  ,r+1);
							hashing.set_data(4,i,j,k,p,q+1,r); hashing.set_data(5,i,j,k,p+1,q+1,r); hashing.set_data(6,i,j,k,p+1,q+1,r+1); hashing.set_data(7,i,j,k,p,q+1,r+1);

							tmp.clear(); bool flag = false;
							contouring(tmp, hashing, _phi, org, h, bdry_xm, bdry_xp, bdry_ym, bdry_yp, bdry_zm, bdry_zp, flag);

							for(size_t m=0;m<tmp.size();++m) {

								mesh.v2f.push_back(tmp[m].index[0]);
								mesh.v2f.push_back(tmp[m].index[1]);
								mesh.v2f.push_back(tmp[m].index[2]);
							}	

						}								
					}
				}
			}			
		}
	}	

	if(hashing.numVertices != int(hashing.posVerts.size()) || hashing.numVertices != int(hashing.vert_types.size()))
		throwError("Error in Marching Interface SRL _ Num Vertices is Different");

	mesh.posVerts = hashing.posVerts;
	types = hashing.vert_types;

	// x-wall
	int num = int(mesh.posVerts.size());
	for(size_t i=0;i<x_wall.size();++i) { mesh.v2f.push_back(x_wall[i].index[0]+num); mesh.v2f.push_back(x_wall[i].index[1]+num); mesh.v2f.push_back(x_wall[i].index[2]+num); }

	vectorLibrary::append(mesh.posVerts, x_hashing.posVerts);
	vectorLibrary::append(types, x_hashing.vert_types);

	// y-wall
	num = int(mesh.posVerts.size());
	for(size_t i=0;i<y_wall.size();++i) { mesh.v2f.push_back(y_wall[i].index[0]+num); mesh.v2f.push_back(y_wall[i].index[1]+num); mesh.v2f.push_back(y_wall[i].index[2]+num); }

	vectorLibrary::append(mesh.posVerts, y_hashing.posVerts);
	vectorLibrary::append(types, y_hashing.vert_types);

	// z-wall
	num = int(mesh.posVerts.size());
	for(size_t i=0;i<z_wall.size();++i) { mesh.v2f.push_back(z_wall[i].index[0]+num); mesh.v2f.push_back(z_wall[i].index[1]+num); mesh.v2f.push_back(z_wall[i].index[2]+num); }

	vectorLibrary::append(mesh.posVerts, z_hashing.posVerts);
	vectorLibrary::append(types, z_hashing.vert_types);

	mesh.numVertices = unsigned int(mesh.posVerts.size());
	mesh.numFaces = unsigned int(mesh.v2f.size() / 3);

	log_info << "SRL3D Marching Cube is Done" << endl;
	log_info << "Number of Vertices : " << mesh.numVertices << endl;
	log_info << "Number of Faces : " << mesh.numFaces << endl;

	std::vector<SRLNode3D<Vector3<T>>> normals_subgrid;
	getNormals_SubGrid(normals_subgrid);

	Grid3<Vector3<T>,T> normals_coarse; normals_coarse.resize(phi);
	MacGridLibrary<T>::getNormals(normals_coarse, phi);

	mesh.normals.resize(mesh.posVerts.size());

#ifdef USE_OPENMP	
	int my_end = int(mesh.posVerts.size());
	int i;

#pragma omp parallel for
#else	
	size_t my_end = mesh.posVerts.size();
	size_t i;
#endif
	for(i=0;i<my_end;++i) {	
//		mesh.normals[i] = Vector3<real>(0,0,1); continue;

		if(types[i] == Triangle3D<real>::BDRY_XM)		{ mesh.normals[i] = Vector3<real>(-1,0,0); continue; }
		else if(types[i] == Triangle3D<real>::BDRY_XP)	{ mesh.normals[i] = Vector3<real>(+1,0,0); continue; }
		else if(types[i] == Triangle3D<real>::BDRY_YM)	{ mesh.normals[i] = Vector3<real>(0,-1,0); continue; }
		else if(types[i] == Triangle3D<real>::BDRY_YP)	{ mesh.normals[i] = Vector3<real>(0,+1,0); continue; }
		else if(types[i] == Triangle3D<real>::BDRY_ZM)	{ mesh.normals[i] = Vector3<real>(0,0,-1); continue; }
		else if(types[i] == Triangle3D<real>::BDRY_ZP)	{ mesh.normals[i] = Vector3<real>(0,0,+1); continue; }

		int p, q, r;
		Vector3<real> x;
		SRL<T>::SRL_Lerp lerp;

		get_idx_normalized_coord(p, q, r, x, mesh.posVerts[i]);

		if(indices(p,q,r) != NOT_REFINE) {
			mesh.normals[i] = unitVector(lerp(*this, normals_subgrid[indices(p,q,r)], x));			

			if(::abs(mesh.normals[i].sqrmag() - 1) > 1.0e-2 || _isnan(mesh.normals[i][0]) || _isnan(mesh.normals[i][1]) || _isnan(mesh.normals[i][2])) {
//				std::cout << "Correction of Normal 1 in SRL3D::marchingInterface()" << std::endl;

				mesh.normals[i] = unitVector(normals_coarse(mesh.posVerts[i]));

				if(::abs(mesh.normals[i].sqrmag() - 1) > 1.0e-2 || _isnan(mesh.normals[i][0]) || _isnan(mesh.normals[i][1]) || _isnan(mesh.normals[i][2])) {
//					std::cout << "Correction of Normal 2 in SRL3D::marchingInterface()" << std::endl;

					mesh.normals[i] = Vector3<real>(0,0,1);
				}
			}			
		}
		else {
			mesh.normals[i] = unitVector(normals_coarse(mesh.posVerts[i]));

			if(::abs(mesh.normals[i].sqrmag() - 1) > 1.0e-2 || _isnan(mesh.normals[i][0]) || _isnan(mesh.normals[i][1]) || _isnan(mesh.normals[i][2])) {
//				std::cout << "Wiered situation in SRL3D::marchingInterface()" << std::endl;
				mesh.normals[i] = Vector3<real>(0,0,1);
			}
		}
	}


//	mesh.init();	
//
//	vector<Triangle3D<T>::TYPE> types;
//
//	vector<Triangle3D<T>>	tmp;
//	Vector3<T>				org, h;
//
//	T _phi[8];
//	bool bdry_xm, bdry_xp, bdry_ym, bdry_yp, bdry_zm, bdry_zp;
//
//	for(index_type i=0;i<phi.nX-1;++i)	{
//		for(index_type j=0;j<phi.nY-1;++j)	{
//			for(index_type k=0;k<phi.nZ-1;++k)	{
//
//				if(isbdry(i,j,k  ) && isbdry(i+1,j,k  ) && isbdry(i,j+1,k  ) && isbdry(i+1,j+1,k  )
//				&& isbdry(i,j,k+1) && isbdry(i+1,j,k+1) && isbdry(i,j+1,k+1) && isbdry(i+1,j+1,k+1) )
//					continue;
//
//				if(indices(i,j,k) == NOT_REFINE) {
//					if(!(phi(i,j,k  ) >  0 && phi(i+1,j,k  ) >  0 && phi(i,j+1,k  ) >  0 && phi(i+1,j+1,k  ) >  0 &&
//						 phi(i,j,k+1) >  0 && phi(i+1,j,k+1) >  0 && phi(i,j+1,k+1) >  0 && phi(i+1,j+1,k+1) >  0 ) 
//					&& !(phi(i,j,k  ) <= 0 && phi(i+1,j,k  ) <= 0 && phi(i,j+1,k  ) <= 0 && phi(i+1,j+1,k  ) <= 0 &&
//						 phi(i,j,k+1) <= 0 && phi(i+1,j,k+1) <= 0 && phi(i,j+1,k+1) <= 0 && phi(i+1,j+1,k+1) <= 0 ) ) {
//
//						cout << i << " " << j << " " << k << endl;
//						cout << isbdry(i,j,k  ) << " " << isbdry(i+1,j,k  ) << " " << isbdry(i,j+1,k  ) << " " << isbdry(i+1,j+1,k  ) << endl;
//						cout << isbdry(i,j,k+1) << " " << isbdry(i+1,j,k+1) << " " << isbdry(i,j+1,k+1) << " " << isbdry(i+1,j+1,k+1) << endl;
//					    cout << phi(i,j,k  ) << " " << phi(i+1,j,k  ) << " " << phi(i,j+1,k  ) << " " << phi(i+1,j+1,k  ) << endl;
//						cout << phi(i,j,k+1) << " " << phi(i+1,j,k+1) << " " << phi(i,j+1,k+1) << " " << phi(i+1,j+1,k+1) << endl;
//						throwError("Error in SRL3D<T,N>::marchingInterface()");
//
//					}					
//
//					if(phi(i,j,k) > 0 || (!isbdry(i,j,k  ) && !isbdry(i+1,j,k  ) && !isbdry(i,j+1,k  ) && !isbdry(i+1,j+1,k  )
//									   && !isbdry(i,j,k+1) && !isbdry(i+1,j,k+1) && !isbdry(i,j+1,k+1) && !isbdry(i+1,j+1,k+1) )
//					  )
//						continue;
//
//					_phi[0] = _phi[1] = _phi[2] = _phi[3] = _phi[4] = _phi[5] = _phi[6] = _phi[7] = -1;
//					h[0] = phi.dx*0.5f; h[1] = phi.dy*0.5f; h[2] = phi.dz*0.5f;
//
//					for(index_type p=0;p<2;++p) {
//						for(index_type q=0;q<2;++q) {
//							for(index_type r=0;r<2;++r) {
//								if(isbdry(i+p,j+q,k+r))
//									continue;
//
//								bdry_xm = bdry_xp = bdry_ym = bdry_yp = bdry_zm = bdry_zp = false;
//
//								if(p==0)	bdry_xp = isbdry(i+1,j+q,k+r);
//								else 		bdry_xm = isbdry(i  ,j+q,k+r);
//
//								if(q==0)	bdry_yp = isbdry(i+p,j+1,k+r);
//								else		bdry_ym = isbdry(i+p,j  ,k+r);
//
//								if(r==0)	bdry_zp = isbdry(i+p,j+q,k+1);
//								else		bdry_zm = isbdry(i+p,j+q,k  );
//
//								org = phi.pos(i,j,k) + Vector3<T>(h[0]*p, h[1]*q, h[2]*r);
//								contouring(tmp, _phi, org, h, bdry_xm, bdry_xp, bdry_ym, bdry_yp, bdry_zm, bdry_zp);
//
//								for(size_t m=0;m<tmp.size();++m) {
//
//									mesh.posVerts.push_back(tmp[m].v[0]);	mesh.v2f.push_back(unsigned int(mesh.posVerts.size()-1)); types.push_back(tmp[m].type);
//									mesh.posVerts.push_back(tmp[m].v[1]);	mesh.v2f.push_back(unsigned int(mesh.posVerts.size()-1)); types.push_back(tmp[m].type);	
//									mesh.posVerts.push_back(tmp[m].v[2]);	mesh.v2f.push_back(unsigned int(mesh.posVerts.size()-1)); types.push_back(tmp[m].type);
//
//								}		
//							}							
//						}
//					}				
//					continue;
//				}
//
//				const int half_N = N/2;
//
//				const SRLNode3D<T,N>& node = nodes[indices(i,j,k)];
//
//				for(int p=0;p<N-1;++p) {
//					for(int q=0;q<N-1;++q) {
//						for(int r=0;r<N-1;++r) {
//
//							if		( p <  half_N && q <  half_N && r <  half_N && isbdry(i  ,j  ,k  ) ) continue;
//							else if	( p >= half_N && q <  half_N && r <  half_N && isbdry(i+1,j  ,k  ) ) continue;
//							else if	( p <  half_N && q >= half_N && r <  half_N && isbdry(i  ,j+1,k  ) ) continue;
//							else if	( p >= half_N && q >= half_N && r <  half_N && isbdry(i+1,j+1,k  ) ) continue;
//
//							else if	( p <  half_N && q <  half_N && r >= half_N && isbdry(i  ,j  ,k+1) ) continue;
//							else if	( p >= half_N && q <  half_N && r >= half_N && isbdry(i+1,j  ,k+1) ) continue;
//							else if	( p <  half_N && q >= half_N && r >= half_N && isbdry(i  ,j+1,k+1) ) continue;
//							else if	( p >= half_N && q >= half_N && r >= half_N && isbdry(i+1,j+1,k+1) ) continue;
//
//							int i_ = (p < half_N) ? 0 : 1;
//							int j_ = (q < half_N) ? 0 : 1;
//							int k_ = (r < half_N) ? 0 : 1;
//
//							bdry_xm = bdry_xp = bdry_ym = bdry_yp = bdry_zm = bdry_zp = false;
//
//							if(p==half_N-1)		bdry_xp = isbdry(i+1,j+j_,k+k_);
//							else if(p==half_N)	bdry_xm = isbdry(i  ,j+j_,k+k_);
//
//							if(q==half_N-1)		bdry_yp = isbdry(i+i_,j+1,k+k_);
//							else if(q==half_N)	bdry_ym = isbdry(i+i_,j  ,k+k_);
//
//							if(r==half_N-1)		bdry_zp = isbdry(i+i_,j+j_,k+1);
//							else if(r==half_N)	bdry_zm = isbdry(i+i_,j+j_,k  );
//
//							org = pos(node, p, q, r);
//							h[0] = h_sub[p]; h[1] = h_sub[q]; h[2] = h_sub[r];
//
//							_phi[0] = node(p,q  ,r);	_phi[1] = node(p+1,q  ,r);	_phi[2] = node(p+1,q  ,r+1);	_phi[3] = node(p,q  ,r+1);
//							_phi[4] = node(p,q+1,r);	_phi[5] = node(p+1,q+1,r);	_phi[6] = node(p+1,q+1,r+1);	_phi[7] = node(p,q+1,r+1);	
//
//							contouring(tmp, _phi, org, h, bdry_xm, bdry_xp, bdry_ym, bdry_yp, bdry_zm, bdry_zp);
//
//							for(size_t m=0;m<tmp.size();++m) {
//
//								mesh.posVerts.push_back(tmp[m].v[0]); mesh.v2f.push_back(unsigned int(mesh.posVerts.size()-1)); types.push_back(tmp[m].type);
//								mesh.posVerts.push_back(tmp[m].v[1]); mesh.v2f.push_back(unsigned int(mesh.posVerts.size()-1)); types.push_back(tmp[m].type);
//								mesh.posVerts.push_back(tmp[m].v[2]); mesh.v2f.push_back(unsigned int(mesh.posVerts.size()-1)); types.push_back(tmp[m].type);
//
//							}	
//
//						}								
//					}
//				}
//			}			
//		}
//	}	
//
//	mesh.numVertices = unsigned int(mesh.posVerts.size());
//	mesh.numFaces = unsigned int(mesh.v2f.size() / 3);
//	
//	std::vector<SRLNode3D<Vector3<T>>> normals_subgrid;
//	getNormals_SubGrid(normals_subgrid);
//
//	Grid3<Vector3<T>,T> normals_coarse; normals_coarse.resize(phi);
//	MacGridLibrary<T>::getNormals(normals_coarse, phi);
//
//	mesh.normals.resize(mesh.posVerts.size());
//	
//#ifdef USE_OPENMP	
//	int my_end = int(mesh.posVerts.size());
//	int i;
//
//#pragma omp parallel for
//#else	
//	size_t my_end = mesh.posVerts.size();
//	size_t i;
//#endif
//	for(i=0;i<my_end;++i) {	
//		if(types[i] == Triangle3D<real>::BDRY_XM)		{ mesh.normals[i] = Vector3<real>(-1,0,0); continue; }
//		else if(types[i] == Triangle3D<real>::BDRY_XP)	{ mesh.normals[i] = Vector3<real>(+1,0,0); continue; }
//		else if(types[i] == Triangle3D<real>::BDRY_YM)	{ mesh.normals[i] = Vector3<real>(0,-1,0); continue; }
//		else if(types[i] == Triangle3D<real>::BDRY_YP)	{ mesh.normals[i] = Vector3<real>(0,+1,0); continue; }
//		else if(types[i] == Triangle3D<real>::BDRY_ZM)	{ mesh.normals[i] = Vector3<real>(0,0,-1); continue; }
//		else if(types[i] == Triangle3D<real>::BDRY_ZP)	{ mesh.normals[i] = Vector3<real>(0,0,+1); continue; }
//
//		int p, q, r;
//		Vector3<real> x;
//		SRL<T>::SRL_Lerp lerp;
//
//		get_idx_normalized_coord(p, q, r, x, mesh.posVerts[i]);
//
//		if(indices(p,q,r) != NOT_REFINE)	mesh.normals[i] = unitVector(lerp(*this, normals_subgrid[indices(p,q,r)], x));			
//		else								mesh.normals[i] = unitVector(normals_coarse(mesh.posVerts[i]));
//	}
}

template<typename T, int N>
void SRL3D<T,N>::glDrawSurfaceMesh(const Material& mat) const {
	mesh.glDraw_SmoothSurface(mat);
}

template<typename T, int N>
void SRL3D<T,N>::glDrawLineStrip(const Material& mat) const {
	mesh.glDrawLineStrip(mat);
}

template<typename T, int N>
void SRL3D<T,N>::glDrawPhi(const RGBA& cp, const RGBA &cn, T band, GLfloat ptsize) const {

	glPointSize(ptsize);

	glBegin(GL_POINTS);
	for(index_type i=1;i<phi.nX-1;i++)	{
		for(index_type j=1;j<phi.nY-1;j++)	{
			for(index_type k=1;k<phi.nZ-1;k++)	{
				T scale = (band - abs(phi(i,j,k)) * hInv) / band;

				if(scale < 0)				continue;
				else if(phi(i,j,k) > 0)		glColor4f(cp[0] * scale, cp[1] * scale, cp[2] * scale, cp[3]);
				else						glColor4f(cn[0] * scale, cn[1] * scale, cn[2] * scale, cn[3]);

				Vector3<T> pos = phi.pos(i,j,k);
				glVertex3f(pos[0], pos[1], pos[2]);
			}			
		}
	}	
	glEnd();
}

template<typename T, int N>
void SRL3D<T,N>::glDrawPhi_SubGrid(const RGBA& cp, const RGBA &cn, T band, GLfloat ptsize) const {

	glPointSize(ptsize);

	Vector3<T> pos0, pos;

	glBegin(GL_POINTS);
	for(size_t i=0;i<nodes.size();++i) {
		pos0 = phi.pos(nodes[i].index[0], nodes[i].index[1], nodes[i].index[2]);
		
		for(index_type p=0;p<N;++p) {
			for(index_type q=0;q<N;++q) {
				for(index_type r=0;r<N;++r) {
					T val = nodes[i].val[p][q][r];
					T scale = (band - abs(val) * hInv) / band;

					if(scale < 0)		continue;
					else if(val > 0)	glColor4f(cp[0] * scale, cp[1] * scale, cp[2] * scale, cp[3]);
					else				glColor4f(cn[0] * scale, cn[1] * scale, cn[2] * scale, cn[3]);

					pos = pos0 + Vector3<T>(pos_sub[p], pos_sub[q], pos_sub[r]);
					glVertex3f(pos[0], pos[1], pos[2]);
				}				
			}			
		}				
	}
	glEnd();
}

template<typename T, int N>
void SRL3D<T,N>::glDrawGrid(const RGBA& c, GLfloat width) const {	
}

template<typename T, int N>
void SRL3D<T,N>::glDrawSubGrid(const RGBA& c, GLfloat ptsize) const {

	glPointSize(ptsize);
	glColor4f( c[0], c[1], c[2], c[3] );

	Vector3<T> pos0, pos;

	glBegin(GL_POINTS);
	for(size_t i=0;i<nodes.size();++i) {
		pos0 = phi.pos(nodes[i].index[0], nodes[i].index[1], nodes[i].index[2]);
		
		for(index_type p=0;p<N;++p) {
			for(index_type q=0;q<N;++q) {
				for(index_type r=0;r<N;++r) {
					pos = pos0 + Vector3<T>(pos_sub[p], pos_sub[q], pos_sub[r]);
					glVertex3f(pos[0], pos[1], pos[2]);
				}				
			}			
		}				
	}
	glEnd();
}

};		// geometryspace