#pragma once

#include "SRL.h"
#include "SRLNode.h"
#include "SRL_Interpolation.h"

namespace geometryspace {

#ifdef _DEBUG
#define __ASSERT(c) assert(c)
#else
#define __ASSERT(c) {}
#endif

template<typename T, int N=SRL_DEFAULT_QUAD_NODE>
class SRL2D : public SRL<T,N>, public MacGrid2D<T> {
public :
	typedef typename SRL<T,N>::index_type index_type;

	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	SRL2D() {}
	SRL2D(const SRL2D& s)	{ set(s);		}
	~SRL2D()				{ release();	}

	void release();

	//////////////////////////////////////////////////////////////////////////
	// Initialization
	void set(index_type nX, index_type nY, T h, T narrowband = 4, bool uniform_quadrature = false);
	void setPhi(const DistanceField2D<T>&);
	void addSource(const DistanceField2D<T>&);

	const SRL2D& operator=(const SRL2D& s) { set(s); }	
	void set(const SRL2D&);
	void set_value(const SRL2D&);

	//////////////////////////////////////////////////////////////////////////
	// Structure 
	void init_subgrid();
	void refine_subgrid();																					// phi needs to be updated properly 
	template<typename D> void refine_subgrid(const Grid2<D,T>&, const D tag_narrow, const D tag_far);		// phi needs to be updated properly 

	template<typename D> static void extract_narrowband(Grid2<D,T>& tagGrid, const Grid2<T,T>& phi, const T band, const D tag_narrow, const D tag_far);

	void modify_coarse_phi_except_near_interface();

	void synchronize_grid_values();	
	void make_node_pointer();

	const SRLNode2D<T,N>* neighbor(const SRLNode2D<T,N>&, int x, int y) const;			// x, y = -1, 0, 1
	const SRLNode2D<T,N>* neighbor(const SRLNode2D<T,N>&, const int x[2]) const;		// x[0], x[1] = -1, 0, 1

	// Used for getting Derivatives _ 1st upwind, 2nd upwind, ENO, WENO
	// Used for Re-initialization
	// RIGHT = 1 or 0
	template<int DIM, int DIR_X, int DIR_Y, int RIGHT> const T neighbor_coarse_level_1(const SRLNode2D<value_type>&, const unsigned int p[2]) const;
	template<int DIM, int DIR_X, int DIR_Y, int RIGHT> const T neighbor_coarse_level_2(const SRLNode2D<value_type>&, const unsigned int p[2]) const;

	template<int DIM, int BDRY, int DIR_X, int DIR_Y, int RIGHT> void neighbor_1(const SRLNode2D<value_type>&, const unsigned int p[2], T v[1], T h[1]) const;
	template<int DIM, int BDRY, int BDRY_INNER, int DIR_X, int DIR_Y, int RIGHT> const bool neighbor_2(const SRLNode2D<value_type>&, const unsigned int p[2], T v[2], T h[2]) const;

	template<int DIM, int DIR_X, int DIR_Y, int RIGHT> const bool neighbor_3(const SRLNode2D<value_type>&, const unsigned int p[2], T v[3]) const;

	void neighbor_x_right_1(const SRLNode2D<value_type>& n, const unsigned int p[2], T v[1], T h[1]) const { neighbor_1<0,N-1,+1, 0, 1>(n, p, v, h); }
	void neighbor_x_left_1(const SRLNode2D<value_type>& n, const unsigned int p[2], T v[1], T h[1]) const  { neighbor_1<0,  0,-1, 0, 0>(n, p, v, h); }
	void neighbor_y_right_1(const SRLNode2D<value_type>& n, const unsigned int p[2], T v[1], T h[1]) const { neighbor_1<1,N-1, 0,+1, 1>(n, p, v, h); }
	void neighbor_y_left_1(const SRLNode2D<value_type>& n, const unsigned int p[2], T v[1], T h[1]) const  { neighbor_1<1,  0, 0,-1, 0>(n, p, v, h); }

	void neighbor_x_right_2(const SRLNode2D<value_type>& n, const unsigned int p[2], T v[2], T h[2]) const { neighbor_2<0,N-1,N-2,+1, 0, 1>(n, p, v, h); }
	void neighbor_x_left_2(const SRLNode2D<value_type>& n, const unsigned int p[2], T v[2], T h[2]) const  { neighbor_2<0,  0,  1,-1, 0, 0>(n, p, v, h); std::swap(v[0], v[1]); std::swap(h[0], h[1]); }
	void neighbor_y_right_2(const SRLNode2D<value_type>& n, const unsigned int p[2], T v[2], T h[2]) const { neighbor_2<1,N-1,N-2, 0,+1, 1>(n, p, v, h); }
	void neighbor_y_left_2(const SRLNode2D<value_type>& n, const unsigned int p[2], T v[2], T h[2]) const  { neighbor_2<1,  0,  1, 0,-1, 0>(n, p, v, h); std::swap(v[0], v[1]); std::swap(h[0], h[1]); }

	const bool neighbor_x_right_2(const SRLNode2D<value_type>& n, const unsigned int p[2], T v[2]) const { 
		T h[2];
		return neighbor_2<0,N-1,N-2,+1, 0, 1>(n, p, v, h); 
	}
	const bool neighbor_x_left_2(const SRLNode2D<value_type>& n, const unsigned int p[2], T v[2]) const  { 
		T h[2];
		const bool flag = neighbor_2<0,  0,  1,-1, 0, 0>(n, p, v, h); 
		std::swap(v[0], v[1]);  
		return flag;
	}
	const bool neighbor_y_right_2(const SRLNode2D<value_type>& n, const unsigned int p[2], T v[2]) const { 
		T h[2];
		return neighbor_2<1,N-1,N-2, 0,+1, 1>(n, p, v, h); 
	}
	const bool neighbor_y_left_2(const SRLNode2D<value_type>& n, const unsigned int p[2], T v[2]) const  { 
		T h[2];
		const bool flag = neighbor_2<1,  0,  1, 0,-1, 0>(n, p, v, h); 
		std::swap(v[0], v[1]); 
		return flag;
	}

	const bool neighbor_x_right_3(const SRLNode2D<value_type>& n, const unsigned int p[2], T v[3]) const { return neighbor_3<0,1,0,1>(n,p,v); }
	const bool neighbor_x_left_3(const SRLNode2D<value_type>& n, const unsigned int p[2], T v[3]) const { 
		if(neighbor_3<0,-1,0,0>(n,p,v))	{ std::swap(v[0],v[2]); return true; }
		else							return false;
	}
	const bool neighbor_y_right_3(const SRLNode2D<value_type>& n, const unsigned int p[2], T v[3]) const { return neighbor_3<1,0,1,1>(n,p,v); }
	const bool neighbor_y_left_3(const SRLNode2D<value_type>& n, const unsigned int p[2], T v[3]) const {
		if(neighbor_3<1,0,-1,0>(n,p,v))	{ std::swap(v[0],v[2]); return true; }
		else							return false;
	}

	//////////////////////////////////////////////////////////////////////////
	// Function for Interpolation
	const bool is_exist_node(index_type i, index_type j) const { return indices(i,j) != NOT_REFINE; }

	const SRLNode2D<T,N>&	node(index_type i, index_type j) const	{ return nodes[indices(i,j)]; }
	SRLNode2D<T,N>&			node(index_type i, index_type j)		{ return nodes[indices(i,j)]; }

	const Vector2<T> pos(const SRLNode2D<T,N>&, index_type p, index_type q) const;

	void get_local_values(T v[2][2], index_type i, index_type j, const Vector2<T>& x) const;					// x : normalized coordinate
	void get_idx_normalized_coord(index_type& i, index_type& j, Vector2<T>& x, const Vector2<T>& pos) const;	// x : normalized coordinate
	
	void get_updating_subgrid_index(Vector2<unsigned int>&, Vector2<unsigned int>&, const SRLNode2D<T,N>&) const;
	void get_updating_subgrid_index(Vector2<unsigned int>&, Vector2<unsigned int>&, const index_type i[2]) const;	

	//////////////////////////////////////////////////////////////////////////
	// Utility
	void createRenderingPrimitive();
	void createRenderingPrimitive(const Grid2<bool,T>& isbdry);

	void marchingTriangle(vector<Triangle2D<T>>&, const Grid2<bool,T>&);
	void marchingInterface(vector<Segment2D<T>>&, const Grid2<bool,T>&);

	void glDrawLiquid(const RGBA&) const;
	void glDrawInterfaceEdges(const RGBA &color = RGBA::blue(), GLfloat width = 3.0f) const;

	void glDrawPhi(const RGBA& cp, const RGBA &cn, T band, GLfloat ptsize) const;	
	void glDrawPhi_SubGrid(const RGBA& cp, const RGBA &cn, T band, GLfloat ptsize) const;
	void glDrawGrid(const RGBA&, GLfloat width) const;
	void glDrawSubGrid(const RGBA&, GLfloat ptsize) const;

	const T getTotalVolume() const;
	const T getTotalVolume(const Grid2<bool,T>& isbdry) const;

	void print_debug_info() const;

	virtual const bool exportBinary(std::ofstream&) const;
	virtual const bool importBinary(std::ifstream&);

	virtual const bool exportBinary_Compressed(std::ofstream&) const;
	virtual const bool importBinary_Compressed(std::ifstream&);

public : 
	Grid2<T,T>							phi;
	Grid2<unsigned int, T>				indices;		// refer index of node in (vector<SRLNode2D> nodes)
														// if not exist, -1

	std::vector<SRLNode2D<T,N>>			nodes;

	// Interpolation
	typename SRL::SRL_Coarse_Lerp		c_lerp;

	//////////////////////////////////////////////////////////////////////////
	// Variable for Visualization	
	mutable std::vector<Triangle2D<T>>	renderingTriangles;
	mutable std::vector<Segment2D<T>>	renderingSegments;
};

template<typename T, int N> 
const bool SRL2D<T,N>::exportBinary(std::ofstream& ofs) const {
	SRL<T,N>::exportBinary(ofs);
	MacGrid2D<T>::exportBinary(ofs);

	phi.exportBinary(ofs);
	indices.exportBinary(ofs);

	unsigned int _size = unsigned int(nodes.size());
	::exportBinary(ofs, _size);

	for(unsigned int i=0;i<_size;++i)
		nodes[i].exportBinary(ofs);	

	return true;
}

template<typename T, int N> 
const bool SRL2D<T,N>::importBinary(std::ifstream& ifs) {
	SRL<T,N>::importBinary(ifs);
	MacGrid2D<T>::importBinary(ifs);

	phi.importBinary(ifs);
	indices.importBinary(ifs);

	unsigned int _size;
	::importBinary(ifs, _size);	

	nodes.resize(_size);
	for(unsigned int i=0;i<_size;++i)
		nodes[i].importBinary(ifs);

	return true;
}

template<typename T, int N> 
const bool SRL2D<T,N>::exportBinary_Compressed(std::ofstream& ofs) const {
	SRL<T,N>::exportBinary(ofs);
	MacGrid2D<T>::exportBinary(ofs);

	Grid2<bool, T> is_negative; is_negative.resize(phi);
	for(index_type i=0;i<is_negative.size();++i) {
		if(phi(i) < 0)	is_negative(i) = true;
		else			is_negative(i) = false;
	}
	is_negative.exportBinary(ofs);	

	std::vector<SRLNode2D<T,N>>	saving_nodes;
	for(std::size_t i=0;i<nodes.size();++i) {
		if(nodes[i].has_interface())
			saving_nodes.push_back(nodes[i]);
	}	

	unsigned int _size = unsigned int(saving_nodes.size());
	::exportBinary(ofs, _size);	

	for(unsigned int i=0;i<_size;++i)
		saving_nodes[i].exportBinary(ofs);

	return true;
}

template<typename T, int N> 
const bool SRL2D<T,N>::importBinary_Compressed(std::ifstream& ifs) {
	SRL<T,N>::importBinary(ifs);
	MacGrid2D<T>::importBinary(ifs);

	Grid2<bool, T> is_negative; 
	is_negative.importBinary(ifs);	

	phi.resize(is_negative);
	indices.resize(is_negative);	indices.set(NOT_REFINE);

	for(index_type i=0;i<is_negative.size();++i) {
		if(is_negative(i))	phi(i) = - (narrowbandWidth + 1);
		else				phi(i) = + (narrowbandWidth + 1);
	}	

	unsigned int _size;
	::importBinary(ifs, _size);	

	nodes.resize(_size);
	for(unsigned int i=0;i<_size;++i) {
		nodes[i].importBinary(ifs);
		indices(nodes[i].index[0], nodes[i].index[1]) = i;
	}

	synchronize_grid_values();

	return true;
}

template<typename T, int N> 
void SRL2D<T,N>::release() { phi.release(); indices.release(); nodes.clear(); renderingTriangles.clear(); }

template<typename T, int N> 
void SRL2D<T,N>::set(index_type nX, index_type nY, T h, T narrowband, bool uniform_quadrature) {
	SRL::set(h, narrowband, uniform_quadrature);
	MacGrid2D::set(nX, nY, h);
	
	const Vector2<T> cOffset(T(-0.5), T(-0.5));
	phi.init(nX+2, nY+2, h, h, cOffset);		phi.zero();
	indices.init(nX+2, nY+2, h, h, cOffset);	indices.set(NOT_REFINE);
}

template<typename T, int N> 
void SRL2D<T,N>::set(const SRL2D<T,N>& src) {
	SRL::set(src.h, src.narrowband, src.uniform_quadrature);
	MacGrid2D<T>::set(src.nX, src.nY, src.h);

	phi.resize(src.phi);			phi = src.phi;
	indices.resize(src.indices);	indices = src.indices;
	nodes = src.nodes;

	renderingTriangles.clear();
}

// Need to Parallelize SRL
template<typename T, int N> 
void SRL2D<T,N>::set_value(const SRL2D<T,N>& src) {
	SRL<T>::SRL_Coarse_Lerp	coarse_lerp;

	phi = src.phi;

#ifdef USE_OPENMP	
	int my_end = int(nodes.size());
	int i;

#pragma omp parallel for
#else	
	size_t my_end = nodes.size();
	size_t i;
#endif
	for( i = 0; i < my_end; ++i ) {			
		const index_type ii = nodes[i].index[0], jj = nodes[i].index[1];
		
		if(src.indices(ii,jj) == NOT_REFINE) {
			Vector2<T> x;

			for(int p=0;p<N;++p) {
				for(int q=0;q<N;++q) {
					x[0] = r_sub[p];
					x[1] = r_sub[q];

					nodes[i].val[p][q] = coarse_lerp(src,ii,jj,x);
				}
			}
		}
		else {
			index_type idx = src.indices(ii,jj);

			for(int p=0;p<N;++p) 
				for(int q=0;q<N;++q) 
					nodes[i].val[p][q] = src.nodes[idx].val[p][q];
		}
	}
}
 
template<typename T, int N>
void SRL2D<T,N>::setPhi(const DistanceField2D<T>& field) {
	for(index_type i=0;i<phi.nX;++i)	{
		for(index_type j=0;j<phi.nY;++j)	{
			Vector2<T> pos = phi.pos(i,j);
			phi(i,j) = field.getPhi(pos);
		}
	}

	refine_subgrid();

	Vector2<T> pos0, pos;
	for(size_t i=0;i<nodes.size();++i) {
		pos0 = phi.pos(nodes[i].index[0],nodes[i].index[1]);
		
		for(int qi=0;qi<N;++qi) {
			for(int qj=0;qj<N;++qj) {
				pos[0] = pos0[0] + pos_sub[qi];
				pos[1] = pos0[1] + pos_sub[qj];

				T val = field.getPhi(pos);

				nodes[i](qi, qj) = val;
			}
		}
	}
}

template<typename T, int N>
void SRL2D<T,N>::addSource(const DistanceField2D<T>& field) {
	for(index_type i=0;i<phi.nX;++i)	{
		for(index_type j=0;j<phi.nY;++j)	{
			Vector2<T> pos = phi.pos(i,j);
			phi(i,j) = MIN(phi(i,j), field.getPhi(pos));
		}
	}

	std::vector<index_type> add_idx;

	for( index_type i = 0; i < phi.nX-1; ++i ) {			
		for(index_type j=0;j<phi.nY-1;++j) {
			if(indices(i,j) != NOT_REFINE) 
				continue;

			if(phi(i,j) > 0 && phi(i+1,j) > 0 && phi(i,j+1) > 0 && phi(i+1,j+1) > 0)
				continue;

			add_idx.push_back(i);
			add_idx.push_back(j);
		}
	}	

	std::vector<SRLNode2D<T,N>>	backup;
	backup = nodes;

	nodes.clear();
	nodes.resize(backup.size() + add_idx.size() / 2);

	for(std::size_t i=0;i<backup.size();++i) {
		nodes[i].index[0] = backup[i].index[0];	
		nodes[i].index[1] = backup[i].index[1];	

		for(int p=0;p<N;++p) {
			for(int q=0;q<N;++q) {
				Vector2<T> _pos = pos(nodes[i], p, q);
				nodes[i].val[p][q] = MIN(backup[i].val[p][q], field.getPhi(_pos));
			}
		}
	}

	std::size_t idx = 0;
	for(std::size_t i=backup.size();i<nodes.size();++i) {
		indices(add_idx[idx], add_idx[idx+1]);
		nodes[i].index[0] = add_idx[idx];
		nodes[i].index[1] = add_idx[idx+1]; idx += 2;

		for(int p=0;p<N;++p) {
			for(int q=0;q<N;++q) {
				Vector2<T> _pos = pos(nodes[i], p, q);
				nodes[i].val[p][q] = field.getPhi(_pos);
			}
		}
	}
}

template<typename T, int N>
void SRL2D<T,N>::init_subgrid() {
	indices.set(NOT_REFINE);
	nodes.clear(); 
}

template<typename T, int N>
void SRL2D<T,N>::refine_subgrid() {	
	extract_narrowband(indices, phi, narrowbandWidth, unsigned int(NEED_REFINE), unsigned int(NOT_REFINE));
	refine_subgrid(indices, unsigned int(NEED_REFINE), unsigned int(NOT_REFINE));
}

// Need to Parallelize SRL
template<typename T, int N> template<typename D> 
void SRL2D<T,N>::refine_subgrid(const Grid2<D,T>& tag, const D tag_narrow, const D tag_far) {
	if(&tag != &indices)
		indices.set(NOT_REFINE);

	unsigned int n_refine = 0;
	for(index_type i=0;i<tag.nX;++i) for(index_type j=0;j<tag.nY;++j) { if(tag(i,j) == tag_narrow) { indices(i,j) = n_refine++; } }

	nodes.clear();
	nodes.resize(n_refine);

#ifdef USE_OPENMP	
	int my_end = int(phi.nX);
	int i;

#pragma omp parallel for
#else	
	index_type my_end = phi.nX;
	index_type i;
#endif
	for( i = 0; i < my_end; ++i ) {			
		for(index_type j=0;j<phi.nY;++j)	{
			if(tag(i,j) == tag_far) {
//				phi(i,j) = SIGN(phi(i,j)) * (narrowbandWidth + h);
				continue;
			}

			nodes[indices(i,j)].set(i,j);
		}
	}	
}

template<typename T, int N>
void SRL2D<T,N>::modify_coarse_phi_except_near_interface() {
#ifdef USE_OPENMP	
	int my_end = int(phi.nX-1);
	int i;

#pragma omp parallel for
#else	
	index_type my_end = phi.nX-1;
	index_type i;
#endif
	for( i = 1; i < my_end; ++i ) {			
		int sign;

		for(index_type j=1;j<phi.nY-1;++j)	{
			sign = SIGN(phi(i,j));

			if(sign != SIGN(phi(i-1,j)) || sign != SIGN(phi(i+1,j)) 
			|| sign != SIGN(phi(i,j-1)) || sign != SIGN(phi(i,j+1)))
				continue;

			phi(i,j) = sign * (narrowbandWidth + h);
		}
	}	

	phi.setNeumannBoundary();
}


template<typename T, int N>
void SRL2D<T,N>::make_node_pointer() {
}

template<typename T, int N>
const SRLNode2D<T,N>* SRL2D<T,N>::neighbor(const SRLNode2D<T,N>& node, int x, int y) const {	// x, y = -1, 0, 1 
	assert(::abs(x) <= 1 && ::abs(y) <= 1);
	if(indices(node.index[0] + x, node.index[1] + y) == NOT_REFINE)	return NULL;
	else															return &nodes[indices(node.index[0]+x,node.index[1]+y)];
}

template<typename T, int N>
const SRLNode2D<T,N>* SRL2D<T,N>::neighbor(const SRLNode2D<T,N>& node, const int x[2]) const {	// x, y = -1, 0, 1 
	return neighbor(node, x[0], x[1]);	
}

template<typename T, int N> template<int DIM, int DIR_X, int DIR_Y, int RIGHT> 
const T SRL2D<T,N>::neighbor_coarse_level_1(const SRLNode2D<value_type>& node, const unsigned int p[2]) const {
	static const int DIR[2] = {DIR_X, DIR_Y};
	
	Vector2<T> x(r_sub[p[0]], r_sub[p[1]]); x[DIM] = T(RIGHT);
	index_type neighbor[2]; ADD2(neighbor, node.index, DIR); 

	return c_lerp(*this, neighbor, x); 	
}

template<typename T, int N> template<int DIM, int DIR_X, int DIR_Y, int RIGHT> 
const T SRL2D<T,N>::neighbor_coarse_level_2(const SRLNode2D<value_type>& node, const unsigned int p[2]) const {
	static const int DIR1[2] = { DIR_X, DIR_Y };
	static const int DIR2[2] = { DIR_X + DIR_X, DIR_Y + DIR_Y };

	Vector2<T> x(r_sub[p[0]], r_sub[p[1]]); x[DIM] = T(RIGHT);
	index_type neighbor[2]; ADD2(neighbor, node.index, DIR2); 

	if(phi.is_valid_range_for_interp<DIM>(neighbor[DIM]) == true)	{ return c_lerp(*this, neighbor, x);									}
	else															{ ADD2(neighbor, node.index, DIR1); return c_lerp(*this, neighbor, x);	}
}	

template<typename T, int N> template<int DIM, int BDRY, int DIR_X, int DIR_Y, int RIGHT> 
void SRL2D<T,N>::neighbor_1(const SRLNode2D<value_type>& node, const unsigned int p[2], T _v[1], T _h[1]) const {
	__ASSERT(node.index[0] != phi.nX-1 && node.index[1] != phi.nY-1);

	static const int DIR[2] = {DIR_X, DIR_Y};

	const SRLNode2D<T>* ptr(NULL);	

	switch(p[DIM]) {
		case BDRY :	if((ptr = this->neighbor(node,DIR)) == NULL) {	
						_v[0] = neighbor_coarse_level_1<DIM,DIR_X,DIR_Y,RIGHT>(node, p);	_h[0] = h;	
					}
					else { 		
						unsigned int p_neighbor[2] = {p[0], p[1]};	
						switch(RIGHT) {
							case 1 :	p_neighbor[DIM] = 1;	_v[0] = (*ptr)(p_neighbor);	_h[0] = h_sub[0];	break;
							default :	p_neighbor[DIM] = N-2;	_v[0] = (*ptr)(p_neighbor);	_h[0] = h_sub[N-2];	break;										
						}
					}

					break;				

		default :			
					unsigned int p_neighbor[2] = {p[0], p[1]};	
					
					ADD2(p_neighbor, p, DIR);	_v[0] = node(p_neighbor);						
					switch(RIGHT) {
						case 1 :	_h[0] = h_sub[p[DIM]];			break;
						default :	_h[0] = h_sub[p_neighbor[DIM]];	break;
					}						
	}
}

template<typename T, int N> template<int DIM, int BDRY, int BDRY_INNER, int DIR_X, int DIR_Y, int RIGHT> 
const bool SRL2D<T,N>::neighbor_2(const SRLNode2D<value_type>& node, const unsigned int p[2], T _v[2], T _h[2]) const {
	__ASSERT(node.index[0] != phi.nX-1 && node.index[1] != phi.nY-1);
		
	static const int DIR[2] =	{ DIR_X		  , DIR_Y		};
	static const int DIR2[2] =	{ DIR_X+DIR_X , DIR_Y+DIR_Y	};

	const SRLNode2D<T>* ptr(NULL);
	unsigned int p_neighbor[2] = {p[0], p[1]};

	switch(p[DIM]) {
		case BDRY_INNER :						
			ADD2(p_neighbor, p, DIR);	_v[0] = node(p_neighbor);	
			switch(RIGHT) {
				case 1 :	_h[0] = h_sub[p[DIM]];			break;
				default :	_h[0] = h_sub[p_neighbor[DIM]];	break;										
			}	

			if((ptr = this->neighbor(node,DIR)) == NULL) {	
				_v[1] = neighbor_coarse_level_1<DIM,DIR_X,DIR_Y,RIGHT>(node, p);	_h[1] = h; 		
				return false;
			}
			else {	
				switch(RIGHT) {
					case 1 :	p_neighbor[DIM] = 1;	_v[1] = (*ptr)(p_neighbor);	_h[1] = h_sub[0];	break;
					default :	p_neighbor[DIM] = N-2;	_v[1] = (*ptr)(p_neighbor);	_h[1] = h_sub[N-2];	break;										
				}
				return true;
			}			

			break;
							
		case BDRY		:	
			if((ptr = this->neighbor(node,DIR)) == NULL) {	
				_v[0] = neighbor_coarse_level_1<DIM,DIR_X,DIR_Y,RIGHT>(node, p);	_h[0] = h; 
				_v[1] = neighbor_coarse_level_2<DIM,DIR_X,DIR_Y,RIGHT>(node, p);	_h[1] = h; 
				return false;
			}	
			else {
				switch(RIGHT) {
					case 1 :	
						p_neighbor[DIM] = 1; _v[0] = (*ptr)(p_neighbor); _h[0] = h_sub[0];	
						p_neighbor[DIM] = 2; _v[1] = (*ptr)(p_neighbor); _h[1] = h_sub[1];							
						break;

					default :	
						p_neighbor[DIM] = N-2; _v[0] = (*ptr)(p_neighbor); _h[0] = h_sub[N-2];	
						p_neighbor[DIM] = N-3; _v[1] = (*ptr)(p_neighbor); _h[1] = h_sub[N-3];							
						break;
				}		
				return true;
			}

			break;

		default			:
			ADD2(p_neighbor, p, DIR);	_v[0] = node(p_neighbor);									
			switch(RIGHT) {
				case 1 :	_h[0] = h_sub[p[DIM]];			break;
				default :	_h[0] = h_sub[p_neighbor[DIM]];	break;	
			}	

			ADD2(p_neighbor, p, DIR2);	_v[1] = node(p_neighbor);									
			switch(RIGHT) {
				case 1 :	_h[1] = h_sub[p[DIM]];			break;
				default :	_h[1] = h_sub[p_neighbor[DIM]];	break;	
			}	
			return true;
	}	
}

template<typename T, int N> template<int DIM, int DIR_X, int DIR_Y, int RIGHT> 
const bool SRL2D<T,N>::neighbor_3(const SRLNode2D<value_type>& node, const unsigned int p[2], T _v[3]) const {
	static const int DIR1[2] = { DIR_X    , DIR_Y     };
	static const int DIR2[2] = { DIR_X * 2, DIR_Y * 2 };
	static const int DIR3[2] = { DIR_X * 3, DIR_Y * 3 };

	const SRLNode2D<T>* ptr(NULL);
	unsigned int p_neighbor[2];

	ADD2(p_neighbor, p, DIR3);
	if(0 <= p_neighbor[DIM] && p_neighbor[DIM] < N) {
		ADD2(p_neighbor, p, DIR1); _v[0] = node(p_neighbor); 
		ADD2(p_neighbor, p, DIR2); _v[1] = node(p_neighbor); 
		ADD2(p_neighbor, p, DIR3); _v[2] = node(p_neighbor); 
		return true;
	}
	
	if( (node.index[DIM] == 0 && RIGHT != 1) || (node.index[DIM] == phi.dim[DIM] - 2 && RIGHT == 1) ) {
		ADD2(p_neighbor, p, DIR1); _v[0] = node(p_neighbor); 
		ADD2(p_neighbor, p, DIR2); _v[1] = node(p_neighbor); 
								   _v[2] = _v[1];								   
		return true;
	}

	if((ptr = this->neighbor(node,DIR1)) == NULL) { return false; }

	ADD2(p_neighbor, p, DIR2);
	if(0 <= p_neighbor[DIM] && p_neighbor[DIM] < N) {
		ADD2(p_neighbor, p, DIR1); _v[0] = node(p_neighbor); 
		ADD2(p_neighbor, p, DIR2); _v[1] = node(p_neighbor); 
		switch(RIGHT) {
			case 1 :	p_neighbor[DIM] = 1;	_v[2] = (*ptr)(p_neighbor); break;
			default :	p_neighbor[DIM] = N-2;	_v[2] = (*ptr)(p_neighbor); break;
		}
		return true;
	}
	
	ADD2(p_neighbor, p, DIR1);
	if(0 <= p_neighbor[DIM] && p_neighbor[DIM] < N) {
		ADD2(p_neighbor, p, DIR1); _v[0] = node(p_neighbor);
		switch(RIGHT) {
			case 1 :	
				p_neighbor[DIM] = 1; _v[1] = (*ptr)(p_neighbor); 
				p_neighbor[DIM] = 2; _v[2] = (*ptr)(p_neighbor); 
				break;

			default :	
				p_neighbor[DIM] = N-2; _v[1] = (*ptr)(p_neighbor);
				p_neighbor[DIM] = N-3; _v[2] = (*ptr)(p_neighbor);
				break;
		}	
		return true;
	}
	
	switch(RIGHT) {
		case 1 :	
			p_neighbor[DIM] = 1; _v[0] = (*ptr)(p_neighbor); 
			p_neighbor[DIM] = 2; _v[1] = (*ptr)(p_neighbor); 
			p_neighbor[DIM] = 3; _v[2] = (*ptr)(p_neighbor); 
			break;

		default :	
			p_neighbor[DIM] = N-2; _v[0] = (*ptr)(p_neighbor); 
			p_neighbor[DIM] = N-3; _v[1] = (*ptr)(p_neighbor); 
			p_neighbor[DIM] = N-4; _v[2] = (*ptr)(p_neighbor); 
			break;
	}	
	return true;
}

// Need to Parallelize SRL
template<typename T, int N> template<typename D> 
void SRL2D<T,N>::extract_narrowband(Grid2<D,T>& tagGrid, const Grid2<T,T>& phi, const T band, const D tag_narrow, const D tag_far) {

#ifdef USE_OPENMP	
	int my_end = int(phi.nX-1);
	int i;

#pragma omp parallel for
#else	
	index_type my_end = phi.nX-1;
	index_type i;
#endif
	for( i = 0; i < my_end; ++i ) {			
		for(index_type j=0;j<phi.nY-1;++j)	{
			if( ( -band < phi(i  ,j  ) && phi(i  ,j  ) < band ) 
			||	( -band < phi(i+1,j  ) && phi(i+1,j  ) < band ) 
			||	( -band < phi(i  ,j+1) && phi(i  ,j+1) < band ) 
			||	( -band < phi(i+1,j+1) && phi(i+1,j+1) < band ) ) 
			
			{ tagGrid(i,j) = tag_narrow; }

			else 
			
			{ tagGrid(i,j) = tag_far; }
		}
	}
}

template<typename T, int N>
const Vector2<T> SRL2D<T,N>::pos(const SRLNode2D<T,N>& node, index_type p, index_type q) const {
	return phi.pos(node.index[0],node.index[1]) + Vector2<T>(pos_sub[p], pos_sub[q]);
}

template<typename T, int N>
void SRL2D<T,N>::get_local_values(T v[2][2], index_type i, index_type j, const Vector2<T>& x) const {
	index_type p = get_sub_index(x[0]);
	index_type q = get_sub_index(x[1]);

	v[0][0] = nodes[indices(i,j)](p,q);
	v[1][0] = nodes[indices(i,j)](p+1,q);
	v[0][1] = nodes[indices(i,j)](p,q+1);
	v[1][1] = nodes[indices(i,j)](p+1,q+1);
}

template<typename T, int N>
void SRL2D<T,N>::get_idx_normalized_coord(index_type& i, index_type& j, Vector2<T>& x, const Vector2<T>& pos) const {
	phi.get_idx_normalized_coord(i,j,x[0],x[1],pos);
}

template<typename T, int N>
void SRL2D<T,N>::get_updating_subgrid_index(Vector2<unsigned int>& rs, Vector2<unsigned int>& re, const SRLNode2D<T,N>& node) const {
	get_updating_subgrid_index(rs, re, node.index[0], node.index[1]);
}

template<typename T, int N>
void SRL2D<T,N>::get_updating_subgrid_index(Vector2<unsigned int>& rs, Vector2<unsigned int>& re, const index_type i[2]) const {
	rs[0] = rs[1] = 0;
	re[0] = re[1] = N-1;

	if(indices(i[0]+1,i[1]  ) == NOT_REFINE) re[0] = N;
	if(indices(i[0]  ,i[1]+1) == NOT_REFINE) re[1] = N;
	if(indices(i[0]+1,i[1]+1) == NOT_REFINE) re[0] = re[1] = N;

	if(i[0] == 0)				rs[0] = N / 2;
	else if(i[0] == phi.nX-2)	re[0] = N / 2 + 1;

	if(i[1] == 0)				rs[1] = N / 2;
	else if(i[1] == phi.nY-2)	re[1] = N / 2 + 1;
}

// Need to Parallelize SRL
template<typename T, int N>
void SRL2D<T,N>::synchronize_grid_values() {
//#ifdef USE_OPENMP	
//	int my_end = int(nodes.size());
//	int k;
//
//#pragma omp parallel for
//#else	
	std::size_t my_end = nodes.size();
	std::size_t k;
//#endif
	for( k = 0; k < my_end; ++k ) {			
		const SRLNode2D<T,N>* ptr(NULL);
		SRLNode2D<T,N>& node = nodes[k];		

		index_type i = node.index[0], j = node.index[1];
		if( i != 0 && j != 0 ) 
			phi(i,j) = node.val[0][0];

		if((ptr = neighbor(node,1,0)) != NULL) { for(int q=0;q<N;++q) node.val[N-1][q] = ptr->val[0][q]; } else phi(i+1,j) = node.val[N-1][0];	
		if((ptr = neighbor(node,0,1)) != NULL) { for(int p=0;p<N;++p) node.val[p][N-1] = ptr->val[p][0]; } else phi(i,j+1) = node.val[0][N-1];		
		if((ptr = neighbor(node,1,1)) != NULL) { node.val[N-1][N-1] = ptr->val[0][0];					 } else phi(i+1,j+1) = node.val[N-1][N-1];	

		if( i == 0 )				node.neumann_x_lower();
		else if( i == phi.nX-2 )	node.neumann_x_upper();

		if( j == 0 )				node.neumann_y_lower();
		else if( j == phi.nY-2 )	node.neumann_y_upper();
	}

	phi.setNeumannBoundary();
}

template<typename T, int N>
void SRL2D<T,N>::print_debug_info() const {
	T sum_phi = 0;
	for(index_type i=0;i<phi.size();++i)
		sum_phi += phi(i);

	T sum_refined = 0;
	for(std::size_t i=0;i<nodes.size();++i) 
		for(int p=0;p<N;++p) 
			for(int q=0;q<N;++q) 
				sum_refined += nodes[i].val[p][q];

	char str[256];
	sprintf_s(str, "SRL Coarse Grid Sum _ %30.20g\n", sum_phi);		log_info << str;
	sprintf_s(str, "SRL Fine Grid Sum _ %30.20g\n", sum_refined);	log_info << str;	
}

#undef __ASSERT

};		// namespace geometryspace

#include "SRL2D_Utility.h"