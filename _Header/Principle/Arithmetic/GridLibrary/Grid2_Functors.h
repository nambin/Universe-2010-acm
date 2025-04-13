#pragma once

#include "../../../Basis/Basis/blocked_range.h"
#include "../../../Basis/Basis/blocked_range2d.h"
#include "../../../Basis/Basis/blocked_range3d.h"
#include "../../../Basis/Basis/Environments.h"
#include "GridLibrary.h"

namespace arithmetic {

template<typename T, typename R> class Grid2;

template<typename T, typename R>
void add(Grid2<T,R>& ret, T c0, const Grid2<T,R>& g0, T c1, const Grid2<T,R>& g1) {
	typedef typename Grid2<T,R>::index_type index_type;

#ifdef USE_OPENMP
	int my_end = int(ret.nX);
	int i;

#pragma omp parallel for 
#else
	index_type my_end = ret.nX;
	index_type i;
#endif
	for( i = 0; i < my_end; ++i ) {
		for(index_type j=0;j<ret.nY;++j) {
			ret(i,j) = c0 * g0(i,j) + c1 * g1(i,j);
		}
	}
}

template<typename T, typename R>
class Grid2_XNeumann {
	Grid2<T,R> &grid;

public : 
	typedef typename Grid2<T,R>::index_type index_type;

	Grid2_XNeumann(Grid2<T,R> &_grid) : grid(_grid) {}
	void operator() ( const basis::blocked_range<index_type>& r ) const {	
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
			grid(0,i) = grid(1,i);
			grid(grid.nX-1,i) = grid(grid.nX-2,i);
		}										
	}	
};

template<typename T, typename R>
class Grid2_YNeumann {
	Grid2<T,R> &grid;

public : 
	typedef typename Grid2<T,R>::index_type index_type;

	Grid2_YNeumann(Grid2<T,R> &_grid) : grid(_grid) {}
	void operator() ( const basis::blocked_range<index_type>& r ) const {	
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
			grid(i,0) = grid(i,1);
			grid(i,grid.nY-1) = grid(i,grid.nY-2);
		}										
	}	
};

template<typename T, typename R>
class Grid2_XNeumann_Scaling {
	Grid2<T,R> &grid;
	T factor;

public : 
	typedef typename Grid2<T,R>::index_type index_type;

	Grid2_XNeumann_Scaling(Grid2<T,R> &_grid, T _fac) : grid(_grid), factor(_fac) {}
	void operator() ( const basis::blocked_range<index_type>& r ) const {	
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
			grid(0,i) = grid(1,i) * factor;
			grid(grid.nX-1,i) = grid(grid.nX-2,i) * factor;
		}										
	}	
};

template<typename T, typename R>
class Grid2_YNeumann_Scaling {
	Grid2<T,R> &grid;
	T factor;

public : 
	typedef typename Grid2<T,R>::index_type index_type;

	Grid2_YNeumann_Scaling(Grid2<T,R> &_grid, T _fac) : grid(_grid), factor(_fac) {}
	void operator() ( const basis::blocked_range<index_type>& r ) const {	
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
			grid(i,0) = grid(i,1) * factor;
			grid(i,grid.nY-1) = grid(i,grid.nY-2) * factor;
		}										
	}	
};


template<typename T, typename R>
class Grid2_XDirichlet {
	Grid2<T,R> &grid;
	T val;

public : 
	typedef typename Grid2<T,R>::index_type index_type;

	Grid2_XDirichlet(Grid2<T,R> &_grid, T _val) : grid(_grid), val(_val) {}
	void operator() ( const basis::blocked_range<index_type>& r ) const {	
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
			grid(0,i) = grid(grid.nX-1,i) = val;
		}										
	}	
};

template<typename T, typename R>
class Grid2_YDirichlet {
	Grid2<T,R> &grid;
	T val;

public : 
	typedef typename Grid2<T,R>::index_type index_type;

	Grid2_YDirichlet(Grid2<T,R> &_grid, T _val) : grid(_grid), val(_val) {}
	void operator() ( const basis::blocked_range<index_type>& r ) const {	
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
			grid(i,0) = grid(i,grid.nY-1) = val;
		}										
	}	
};


};		// arithmetic