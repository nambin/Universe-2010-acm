#pragma once

#include "../../../Basis/Basis/Environments.h"

namespace arithmetic {

template<typename T>
class BuildPreconditioner_BlackBlock_Grid2 {
	typedef typename Grid2_Preconditioner_BlockRedBlackICCG<T>::Block Block;

public : 
	const std::vector<Block>& blocks;

	const Grid2<bool,T>&	isBoundary;
	const Grid2x3<T>&		poisson;

	Grid2<T,T>&			preconditioner;	

	BuildPreconditioner_BlackBlock_Grid2
	(const Grid2<bool,T>& _isBdry, const Grid2x3<T>& _poi, Grid2<T,T>& _pre, const std::vector<Block>& _blocks)
	: isBoundary(_isBdry), poisson(_poi), preconditioner(_pre), blocks(_blocks) {}

	void operator()(const basis::blocked_range<grid_index_type>& r) const {
#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for
#else
		grid_index_type my_begin = r.begin();
		grid_index_type my_end = r.end();
		grid_index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			buildPreconditioner(blocks[i]);
		}
	}
	void buildPreconditioner(const Block& block) const {
		T d;
		grid_index_type i, j;

		for(i=block.begin[0];i<block.end[0]-1;++i) {
			for(j=block.begin[1];j<block.end[1]-1;++j) {
				if( !isBoundary(i,j) )	{
					d = poisson(i,j,0) 
						- SQ( poisson(i-1,j,1)*preconditioner(i-1,j) )
						- SQ( poisson(i,j-1,2)*preconditioner(i,j-1) );
					preconditioner(i,j) = 1.0f / sqrtf((float)::abs(d)+0.001f);
				}
			}
		}

		i = block.end[0]-1;
		for(j=block.begin[1];j<block.end[1]-1;++j) {
			if( !isBoundary(i,j) )	{
				d = poisson(i,j,0) 
					- SQ( poisson(i-1,j,1)*preconditioner(i-1,j) )
					- SQ( poisson(i,j-1,2)*preconditioner(i,j-1) )
					- SQ( poisson(i  ,j,1)*preconditioner(i+1,j) );
				preconditioner(i,j) = 1.0f / sqrtf((float)::abs(d)+0.001f);
			}
		}

		j = block.end[1]-1;
		for(i=block.begin[0];i<block.end[0]-1;++i) {
			if( !isBoundary(i,j) )	{
				d = poisson(i,j,0) 
					- SQ( poisson(i-1,j,1)*preconditioner(i-1,j) )
					- SQ( poisson(i,j-1,2)*preconditioner(i,j-1) )
					- SQ( poisson(i  ,j,2)*preconditioner(i,j+1) );
				preconditioner(i,j) = 1.0f / sqrtf((float)::abs(d)+0.001f);
			}
		}

		i = block.end[0]-1;
		j = block.end[1]-1;

		if(!isBoundary(i,j)) {
			d = poisson(i,j,0) 
				- SQ( poisson(i-1,j,1)*preconditioner(i-1,j) )
				- SQ( poisson(i,j-1,2)*preconditioner(i,j-1) )
				- SQ( poisson(i  ,j,1)*preconditioner(i+1,j) )
				- SQ( poisson(i  ,j,2)*preconditioner(i,j+1) );
			preconditioner(i,j) = 1.0f / sqrtf((float)::abs(d)+0.001f);
		}		
	}
};


template<typename T>
class ForwardSubstitution_BlackBlock_Grid2 {
	typedef typename Grid2_Preconditioner_BlockRedBlackICCG<T>::Block Block;

public : 
	const std::vector<Block>& blocks;

	const Grid2<bool,T>&	isBoundary;
	const Grid2x3<T>&		poisson;
	const Grid2<T,T>&	preconditioner;	

	const Grid2<T, T>&	x;
	Grid2<T, T>&		m;

	ForwardSubstitution_BlackBlock_Grid2
	(const Grid2<bool,T>& _isBdry, const Grid2x3<T>& _poi, const Grid2<T,T>& _pre, 
	const std::vector<Block>& _blocks, const Grid2<T,T> &_x, Grid2<T,T> &_m)
	: isBoundary(_isBdry), poisson(_poi), preconditioner(_pre), blocks(_blocks), x(_x), m(_m) {}

	void operator()(const basis::blocked_range<grid_index_type>& r) const {
#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for
#else
		grid_index_type my_begin = r.begin();
		grid_index_type my_end = r.end();
		grid_index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			forwardSubstitution(blocks[i]);
		}
	}
	void forwardSubstitution(const Block& block) const {
		T d;
		grid_index_type i, j;

		for(i=block.begin[0];i<block.end[0]-1;++i) {
			for(j=block.begin[1];j<block.end[1]-1;++j) {
				if( !isBoundary(i,j) )	{
					d = x(i,j) 
						- poisson(i-1,j,1)*preconditioner(i-1,j)*m(i-1,j)
						- poisson(i,j-1,2)*preconditioner(i,j-1)*m(i,j-1);
					m(i,j) = preconditioner(i,j)*d;
				}
			}
		}

		i = block.end[0]-1;
		for(j=block.begin[1];j<block.end[1]-1;++j) {
			if( !isBoundary(i,j) )	{
				d = x(i,j) 
					- poisson(i-1,j,1)*preconditioner(i-1,j)*m(i-1,j)
					- poisson(i,j-1,2)*preconditioner(i,j-1)*m(i,j-1)
					- poisson(i,j  ,1)*preconditioner(i+1,j)*m(i+1,j);
				m(i,j) = preconditioner(i,j)*d;
			}
		}

		j = block.end[1]-1;
		for(i=block.begin[0];i<block.end[0]-1;++i) {
			if( !isBoundary(i,j) )	{
				d = x(i,j) 
					- poisson(i-1,j,1)*preconditioner(i-1,j)*m(i-1,j)
					- poisson(i,j-1,2)*preconditioner(i,j-1)*m(i,j-1)
					- poisson(i,j  ,2)*preconditioner(i,j+1)*m(i,j+1);
				m(i,j) = preconditioner(i,j)*d;
			}
		}

		i = block.end[0]-1;
		j = block.end[1]-1;

		if(!isBoundary(i,j)) {
			d = x(i,j) 
				- poisson(i-1,j,1)*preconditioner(i-1,j)*m(i-1,j)
				- poisson(i,j-1,2)*preconditioner(i,j-1)*m(i,j-1)
				- poisson(i,j  ,1)*preconditioner(i+1,j)*m(i+1,j)
				- poisson(i,j  ,2)*preconditioner(i,j+1)*m(i,j+1);
			m(i,j) = preconditioner(i,j)*d;
		}	
	}
};

template<typename T>
class BackwardSubstitution_BlackBlock_Grid2 {
	typedef typename Grid2_Preconditioner_BlockRedBlackICCG<T>::Block Block;

public : 
	const std::vector<Block>& blocks;

	const Grid2<bool,T>&	isBoundary;
	const Grid2x3<T>&		poisson;
	const Grid2<T,T>&	preconditioner;	

	const Grid2<T, T>&	m;
	Grid2<T, T>&		y;

	BackwardSubstitution_BlackBlock_Grid2
	(const Grid2<bool,T>& _isBdry, const Grid2x3<T>& _poi, const Grid2<T,T>& _pre, 
	const std::vector<Block>& _blocks, const Grid2<T,T> &_m, Grid2<T,T> &_y)
	: isBoundary(_isBdry), poisson(_poi), preconditioner(_pre), blocks(_blocks), m(_m), y(_y) {}

	void operator()(const basis::blocked_range<grid_index_type>& r) const {
#ifdef USE_OPENMP
		int my_begin = int(r.begin());
		int my_end = int(r.end());
		int i;

#pragma omp parallel for
#else
		grid_index_type my_begin = r.begin();
		grid_index_type my_end = r.end();
		grid_index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			backwardSubstitution(blocks[i]);
		}
	}
	void backwardSubstitution(const Block& block) const {
		T d;
		grid_index_type i, j;

		i = block.end[0]-1;
		j = block.end[1]-1;

		if(!isBoundary(i,j)) {
			d = m(i,j);				
			y(i,j) = preconditioner(i,j)*d;
		}	

		i = block.end[0]-1;
		for(j=block.end[1]-2;j>=block.begin[1];--j) {
			if( !isBoundary(i,j) )	{
				d = m(i,j) 					
					- poisson(i,j,2)*preconditioner(i,j)*y(i,j+1);
				y(i,j) = preconditioner(i,j)*d;
			}
		}

		j = block.end[1]-1;
		for(i=block.end[0]-2;i>=block.begin[0];--i) {
			if( !isBoundary(i,j) )	{
				d = m(i,j) 
					- poisson(i,j,1)*preconditioner(i,j)*y(i+1,j);
				y(i,j) = preconditioner(i,j)*d;
			}
		}

		for(i=block.end[0]-2;i>=block.begin[0];--i) {
			for(j=block.end[1]-2;j>=block.begin[1];--j) {
				if( !isBoundary(i,j) )	{
					d = m(i,j) 
						- poisson(i,j,1)*preconditioner(i,j)*y(i+1,j)
						- poisson(i,j,2)*preconditioner(i,j)*y(i,j+1);
					y(i,j) = preconditioner(i,j)*d;
				}
			}
		}
	}
};

};		// namespace arithmetic