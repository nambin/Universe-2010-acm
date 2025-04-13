#pragma once

#include "../../../Basis/Basis/Environments.h"
#include "GridLibrary.h"

namespace arithmetic {

template<typename T>
class BuildPreconditioner_RedBlock_Grid3 {
	typedef typename Grid3_Preconditioner_BlockRedBlackICCG<T>::Block Block;
	typedef typename Grid3_Preconditioner_BlockRedBlackICCG<T>::index_type index_type;

public : 
	const std::vector<Block>&	blocks;

	const Grid3<bool,T>&		isBoundary;
	const Grid3x4<T>&			poisson;

	Grid3<T,T>&					preconditioner;	

	BuildPreconditioner_RedBlock_Grid3
	(const Grid3<bool,T>& _isBdry, const Grid3x4<T>& _poi, Grid3<T,T>& _pre, const std::vector<Block>& _blocks)
	: isBoundary(_isBdry), poisson(_poi), preconditioner(_pre), blocks(_blocks) {}

	void operator()(const basis::blocked_range<index_type>& r) const {
#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for
#else
		index_type my_begin = r.begin();
		index_type my_end = r.end();
		index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			buildPreconditioner(blocks[i]);
		}
	}
	void buildPreconditioner(const Block& block) const {
		T d;
		index_type i,j,k;

		i = block.begin[0];
		j = block.begin[1];
		k = block.begin[2];
		if(!isBoundary(i,j,k)) {
			d = poisson(i,j,k,0);
			preconditioner(i,j,k) = 1.0f / sqrt((float)::abs(d)+0.001f);
		}

		j = block.begin[1];
		k = block.begin[2];
		for(i=block.begin[0]+1;i<block.end[0];++i) {
			if(!isBoundary(i,j,k)) {
				d = poisson(i,j,k,0) 
					- SQ( poisson(i-1,j,k,1)*preconditioner(i-1,j,k) );
				preconditioner(i,j,k) = 1.0f / sqrt((float)::abs(d)+0.001f);
			}
		}

		i = block.begin[0];
		k = block.begin[2];
		for(j=block.begin[1]+1;j<block.end[1];++j) {
			if(!isBoundary(i,j,k)) {
				d = poisson(i,j,k,0) 
					- SQ( poisson(i,j-1,k,2)*preconditioner(i,j-1,k) );
				preconditioner(i,j,k) = 1.0f / sqrt((float)::abs(d)+0.001f);
			}
		}

		i = block.begin[0];
		j = block.begin[1];
		for(k=block.begin[2]+1;k<block.end[2];++k) {
			if(!isBoundary(i,j,k)) {
				d = poisson(i,j,k,0) 
					- SQ( poisson(i,j,k-1,3)*preconditioner(i,j,k-1) );
				preconditioner(i,j,k) = 1.0f / sqrt((float)::abs(d)+0.001f);
			}
		}

		i = block.begin[0];
		for(j=block.begin[1]+1;j<block.end[1];++j) {
			for(k=block.begin[2]+1;k<block.end[2];++k) {
				if(!isBoundary(i,j,k)) {
					d = poisson(i,j,k,0) 
						- SQ( poisson(i,j-1,k,2)*preconditioner(i,j-1,k) )
						- SQ( poisson(i,j,k-1,3)*preconditioner(i,j,k-1) );
					preconditioner(i,j,k) = 1.0f / sqrt((float)::abs(d)+0.001f);
				}
			}
		}

		j = block.begin[1];
		for(i=block.begin[0]+1;i<block.end[0];++i) {
			for(k=block.begin[2]+1;k<block.end[2];++k) {
				if(!isBoundary(i,j,k)) {
					d = poisson(i,j,k,0) 
						- SQ( poisson(i-1,j,k,1)*preconditioner(i-1,j,k) )
						- SQ( poisson(i,j,k-1,3)*preconditioner(i,j,k-1) );
					preconditioner(i,j,k) = 1.0f / sqrt((float)::abs(d)+0.001f);
				}
			}
		}

		k = block.begin[2];
		for(i=block.begin[0]+1;i<block.end[0];++i) {
			for(j=block.begin[1]+1;j<block.end[1];++j) {
				if(!isBoundary(i,j,k)) {
					d = poisson(i,j,k,0) 
						- SQ( poisson(i-1,j,k,1)*preconditioner(i-1,j,k) )
						- SQ( poisson(i,j-1,k,2)*preconditioner(i,j-1,k) );
					preconditioner(i,j,k) = 1.0f / sqrt((float)::abs(d)+0.001f);
				}
			}
		}
	
		for(i=block.begin[0]+1;i<block.end[0];++i) {
			for(j=block.begin[1]+1;j<block.end[1];++j) {
				for(k=block.begin[2]+1;k<block.end[2];++k) {
					if(!isBoundary(i,j,k)) {
						d = poisson(i,j,k,0) 
							- SQ( poisson(i-1,j,k,1)*preconditioner(i-1,j,k) )
							- SQ( poisson(i,j-1,k,2)*preconditioner(i,j-1,k) )
							- SQ( poisson(i,j,k-1,3)*preconditioner(i,j,k-1) );
						preconditioner(i,j,k) = 1.0f / sqrt((float)::abs(d)+0.001f);
					}
				}
			}
		}
	}
};

template<typename T>
class ForwardSubstitution_RedBlock_Grid3 {
	typedef typename Grid3_Preconditioner_BlockRedBlackICCG<T>::Block Block;
	typedef typename Grid3_Preconditioner_BlockRedBlackICCG<T>::index_type index_type;

public : 
	const std::vector<Block>& blocks;

	const bool*				isBoundary;
	const T*				poisson;
	const T*				preconditioner;	

	const T*				x;
	T*						m;

	index_type				nX, nY, nZ;

	ForwardSubstitution_RedBlock_Grid3
	(const Grid3<bool,T>& _isBdry, const Grid3x4<T>& _poi, const Grid3<T,T>& _pre, 
	const std::vector<Block>& _blocks, const Grid3<T,T> &_x, Grid3<T,T> &_m)
	: isBoundary(_isBdry.getptr()), poisson(_poi.getptr()), preconditioner(_pre.getptr()), blocks(_blocks), x(_x.getptr()), m(_m.getptr()),
	nX(_pre.nX), nY(_pre.nY), nZ(_pre.nZ) {}

	void operator()(const basis::blocked_range<index_type>& r) const {
#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for
#else
		index_type my_begin = r.begin();
		index_type my_end = r.end();
		index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			forwardSubstitution(blocks[i]);
		}
	}
	void forwardSubstitution(const Block block) const {
		T d;
		index_type i,j,k;

#define IS(i,j,k) ((k) + (nZ) * ((j) + (nY) * (i)))

		i = block.begin[0];
		j = block.begin[1];
		k = block.begin[2];
		if(!isBoundary[IS(i,j,k)]) {
			d = x[IS(i,j,k)];
			m[IS(i,j,k)] = preconditioner[IS(i,j,k)]*d;
		}

		j = block.begin[1];
		k = block.begin[2];
		for(i=block.begin[0]+1;i<block.end[0];++i) {
			if(!isBoundary[IS(i,j,k)]) {
				d = x[IS(i,j,k)] 
					- poisson[4*IS(i-1,j,k)+1]*preconditioner[IS(i-1,j,k)]*m[IS(i-1,j,k)];
				m[IS(i,j,k)] = preconditioner[IS(i,j,k)]*d;
			}
		}

		i = block.begin[0];
		k = block.begin[2];
		for(j=block.begin[1]+1;j<block.end[1];++j) {
			if(!isBoundary[IS(i,j,k)]) {
				d = x[IS(i,j,k)] 
					- poisson[4*IS(i,j-1,k)+2]*preconditioner[IS(i,j-1,k)]*m[IS(i,j-1,k)];
				m[IS(i,j,k)] = preconditioner[IS(i,j,k)]*d;
			}
		}

		i = block.begin[0];
		j = block.begin[1];
		for(k=block.begin[2]+1;k<block.end[2];++k) {
			if(!isBoundary[IS(i,j,k)]) {
				d = x[IS(i,j,k)] 				
					- poisson[4*IS(i,j,k-1)+3]*preconditioner[IS(i,j,k-1)]*m[IS(i,j,k-1)];
				m[IS(i,j,k)] = preconditioner[IS(i,j,k)]*d;
			}
		}

		i = block.begin[0];
		for(j=block.begin[1]+1;j<block.end[1];++j) {
			for(k=block.begin[2]+1;k<block.end[2];++k) {
				if(!isBoundary[IS(i,j,k)]) {
					d = x[IS(i,j,k)] 					
						- poisson[4*IS(i,j-1,k)+2]*preconditioner[IS(i,j-1,k)]*m[IS(i,j-1,k)]
						- poisson[4*IS(i,j,k-1)+3]*preconditioner[IS(i,j,k-1)]*m[IS(i,j,k-1)];
					m[IS(i,j,k)] = preconditioner[IS(i,j,k)]*d;
				}
			}
		}

		j = block.begin[1];
		for(i=block.begin[0]+1;i<block.end[0];++i) {
			for(k=block.begin[2]+1;k<block.end[2];++k) {
				if(!isBoundary[IS(i,j,k)]) {
					d = x[IS(i,j,k)] 
						- poisson[4*IS(i-1,j,k)+1]*preconditioner[IS(i-1,j,k)]*m[IS(i-1,j,k)]					
						- poisson[4*IS(i,j,k-1)+3]*preconditioner[IS(i,j,k-1)]*m[IS(i,j,k-1)];
					m[IS(i,j,k)] = preconditioner[IS(i,j,k)]*d;
				}
			}
		}

		k = block.begin[2];
		for(i=block.begin[0]+1;i<block.end[0];++i) {
			for(j=block.begin[1]+1;j<block.end[1];++j) {
				if(!isBoundary[IS(i,j,k)]) {
					d = x[IS(i,j,k)] 
						- poisson[4*IS(i-1,j,k)+1]*preconditioner[IS(i-1,j,k)]*m[IS(i-1,j,k)]
						- poisson[4*IS(i,j-1,k)+2]*preconditioner[IS(i,j-1,k)]*m[IS(i,j-1,k)];					
					m[IS(i,j,k)] = preconditioner[IS(i,j,k)]*d;
				}
			}
		}

		for(i=block.begin[0]+1;i<block.end[0];++i) {
			for(j=block.begin[1]+1;j<block.end[1];++j) {
				for(k=block.begin[2]+1;k<block.end[2];++k) {
					if(!isBoundary[IS(i,j,k)]) {
						d = x[IS(i,j,k)] 
							- poisson[4*IS(i-1,j,k)+1]*preconditioner[IS(i-1,j,k)]*m[IS(i-1,j,k)]
							- poisson[4*IS(i,j-1,k)+2]*preconditioner[IS(i,j-1,k)]*m[IS(i,j-1,k)]
							- poisson[4*IS(i,j,k-1)+3]*preconditioner[IS(i,j,k-1)]*m[IS(i,j,k-1)];
						m[IS(i,j,k)] = preconditioner[IS(i,j,k)]*d;
					}
				}
			}
		}
#undef IS
	}
};

template<typename T>
class BackwardSubstitution_RedBlock_Grid3 {
	typedef typename Grid3_Preconditioner_BlockRedBlackICCG<T>::Block Block;
	typedef typename Grid3_Preconditioner_BlockRedBlackICCG<T>::index_type index_type;

public : 
	const std::vector<Block>& blocks;

	const bool*				isBoundary;
	const T*				poisson;
	const T*				preconditioner;	

	const T*				m;
	T*						y;

	index_type				nX, nY, nZ;

	BackwardSubstitution_RedBlock_Grid3
	(const Grid3<bool,T>& _isBdry, const Grid3x4<T>& _poi, const Grid3<T,T>& _pre, 
	const std::vector<Block>& _blocks, const Grid3<T,T> &_m, Grid3<T,T> &_y)
	: isBoundary(_isBdry.getptr()), poisson(_poi.getptr()), preconditioner(_pre.getptr()), blocks(_blocks), m(_m.getptr()), y(_y.getptr()),
	nX(_pre.nX), nY(_pre.nY), nZ(_pre.nZ) {}

	void operator()(const basis::blocked_range<index_type>& r) const {
#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for
#else
		index_type my_begin = r.begin();
		index_type my_end = r.end();
		index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			backwardSubstitution(blocks[i]);
		}
	}
	void backwardSubstitution(const Block block) const {
		T d;
		index_type i, j, k;

#define IS(i,j,k) ((k) + (nZ) * ((j) + (nY) * (i)))

		for(i=block.end[0]-1;i>block.begin[0];--i) {
			for(j=block.end[1]-1;j>block.begin[1];--j) {
				for(k=block.end[2]-1;k>block.begin[2];--k) {
					if(!isBoundary[IS(i,j,k)]) {
						d = m[IS(i,j,k)] 
							- poisson[4*IS(i,j,k)+1]*preconditioner[IS(i,j,k)]*y[IS(i+1,j,k)]	
							- poisson[4*IS(i,j,k)+2]*preconditioner[IS(i,j,k)]*y[IS(i,j+1,k)]
							- poisson[4*IS(i,j,k)+3]*preconditioner[IS(i,j,k)]*y[IS(i,j,k+1)];
						y[IS(i,j,k)] = preconditioner[IS(i,j,k)]*d;
					}
				}				
			}
		}

		i = block.begin[0];
		for(j=block.end[1]-1;j>block.begin[1];--j) {
			for(k=block.end[2]-1;k>block.begin[2];--k) {
				if(!isBoundary[IS(i,j,k)]) {
					d = m[IS(i,j,k)] 
						- poisson[4*IS(i,j,k)+1]*preconditioner[IS(i,j,k)]*y[IS(i+1,j,k)]	
						- poisson[4*IS(i,j,k)+2]*preconditioner[IS(i,j,k)]*y[IS(i,j+1,k)]
						- poisson[4*IS(i,j,k)+3]*preconditioner[IS(i,j,k)]*y[IS(i,j,k+1)]
						- poisson[4*IS(i-1,j,k)+1]*preconditioner[IS(i,j,k)]*y[IS(i-1,j,k)];
					y[IS(i,j,k)] = preconditioner[IS(i,j,k)]*d;
				}
			}				
		}

		j = block.begin[1];
		for(i=block.end[0]-1;i>block.begin[0];--i) {
			for(k=block.end[2]-1;k>block.begin[2];--k) {
				if(!isBoundary[IS(i,j,k)]) {
					d = m[IS(i,j,k)] 
						- poisson[4*IS(i,j,k)+1]*preconditioner[IS(i,j,k)]*y[IS(i+1,j,k)]	
						- poisson[4*IS(i,j,k)+2]*preconditioner[IS(i,j,k)]*y[IS(i,j+1,k)]
						- poisson[4*IS(i,j,k)+3]*preconditioner[IS(i,j,k)]*y[IS(i,j,k+1)]
						- poisson[4*IS(i,j-1,k)+2]*preconditioner[IS(i,j,k)]*y[IS(i,j-1,k)];
					y[IS(i,j,k)] = preconditioner[IS(i,j,k)]*d;
				}
			}				
		}

		k = block.begin[2];
		for(i=block.end[0]-1;i>block.begin[0];--i) {
			for(j=block.end[1]-1;j>block.begin[1];--j) {
				if(!isBoundary[IS(i,j,k)]) {
					d = m[IS(i,j,k)] 
						- poisson[4*IS(i,j,k)+1]*preconditioner[IS(i,j,k)]*y[IS(i+1,j,k)]	
						- poisson[4*IS(i,j,k)+2]*preconditioner[IS(i,j,k)]*y[IS(i,j+1,k)]
						- poisson[4*IS(i,j,k)+3]*preconditioner[IS(i,j,k)]*y[IS(i,j,k+1)]
						- poisson[4*IS(i,j,k-1)+3]*preconditioner[IS(i,j,k)]*y[IS(i,j,k-1)];
					y[IS(i,j,k)] = preconditioner[IS(i,j,k)]*d;
				}
			}				
		}

		i = block.begin[0];
		j = block.begin[1];
		for(k=block.end[2]-1;k>block.begin[2];--k) {
			if(!isBoundary[IS(i,j,k)]) {
				d = m[IS(i,j,k)] 
					- poisson[4*IS(i,j,k)+1]*preconditioner[IS(i,j,k)]*y[IS(i+1,j,k)]	
					- poisson[4*IS(i,j,k)+2]*preconditioner[IS(i,j,k)]*y[IS(i,j+1,k)]
					- poisson[4*IS(i,j,k)+3]*preconditioner[IS(i,j,k)]*y[IS(i,j,k+1)]
					- poisson[4*IS(i-1,j,k)+1]*preconditioner[IS(i,j,k)]*y[IS(i-1,j,k)]
					- poisson[4*IS(i,j-1,k)+2]*preconditioner[IS(i,j,k)]*y[IS(i,j-1,k)];
				y[IS(i,j,k)] = preconditioner[IS(i,j,k)]*d;
			}
		}	

		i = block.begin[0];
		k = block.begin[2];
		for(j=block.end[1]-1;j>block.begin[1];--j) {
			if(!isBoundary[IS(i,j,k)]) {
				d = m[IS(i,j,k)] 
					- poisson[4*IS(i,j,k)+1]*preconditioner[IS(i,j,k)]*y[IS(i+1,j,k)]	
					- poisson[4*IS(i,j,k)+2]*preconditioner[IS(i,j,k)]*y[IS(i,j+1,k)]
					- poisson[4*IS(i,j,k)+3]*preconditioner[IS(i,j,k)]*y[IS(i,j,k+1)]
					- poisson[4*IS(i-1,j,k)+1]*preconditioner[IS(i,j,k)]*y[IS(i-1,j,k)]
					- poisson[4*IS(i,j,k-1)+3]*preconditioner[IS(i,j,k)]*y[IS(i,j,k-1)];
				y[IS(i,j,k)] = preconditioner[IS(i,j,k)]*d;
			}
		}	

		j = block.begin[1];
		k = block.begin[2];
		for(i=block.end[0]-1;i>block.begin[0];--i) {
			if(!isBoundary[IS(i,j,k)]) {
				d = m[IS(i,j,k)] 
					- poisson[4*IS(i,j,k)+1]*preconditioner[IS(i,j,k)]*y[IS(i+1,j,k)]	
					- poisson[4*IS(i,j,k)+2]*preconditioner[IS(i,j,k)]*y[IS(i,j+1,k)]
					- poisson[4*IS(i,j,k)+3]*preconditioner[IS(i,j,k)]*y[IS(i,j,k+1)]
					- poisson[4*IS(i,j-1,k)+2]*preconditioner[IS(i,j,k)]*y[IS(i,j-1,k)]
					- poisson[4*IS(i,j,k-1)+3]*preconditioner[IS(i,j,k)]*y[IS(i,j,k-1)];
				y[IS(i,j,k)] = preconditioner[IS(i,j,k)]*d;
			}
		}	

		i = block.begin[0];
		j = block.begin[1];
		k = block.begin[2];
		if(!isBoundary[IS(i,j,k)]) {
			d = m[IS(i,j,k)] 
				- poisson[4*IS(i,j,k)+1]*preconditioner[IS(i,j,k)]*y[IS(i+1,j,k)]	
				- poisson[4*IS(i,j,k)+2]*preconditioner[IS(i,j,k)]*y[IS(i,j+1,k)]
				- poisson[4*IS(i,j,k)+3]*preconditioner[IS(i,j,k)]*y[IS(i,j,k+1)]
				- poisson[4*IS(i-1,j,k)+1]*preconditioner[IS(i,j,k)]*y[IS(i-1,j,k)]
				- poisson[4*IS(i,j-1,k)+2]*preconditioner[IS(i,j,k)]*y[IS(i,j-1,k)]
				- poisson[4*IS(i,j,k-1)+3]*preconditioner[IS(i,j,k)]*y[IS(i,j,k-1)];
			y[IS(i,j,k)] = preconditioner[IS(i,j,k)]*d;
		}

#undef IS
	}
};

};		// arithmetic