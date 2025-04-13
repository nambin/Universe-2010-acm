#pragma once

#include "../../../Basis/Basis/blocked_range.h"
#include "../../../Basis/Basis/blocked_range2d.h"
#include "../../../Basis/Basis/blocked_range3d.h"
#include "../../../Basis/Basis/Environments.h"
#include "GridLibrary.h"

namespace arithmetic {

template<typename T, typename R> class Grid3;

template<typename T, typename R>
void add(Grid3<T,R>& ret, T c0, const Grid3<T,R>& g0, T c1, const Grid3<T,R>& g1) {
	typedef typename Grid3<T,R>::index_type index_type;

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
			for(index_type k=0;k<ret.nZ;++k) {
				ret(i,j,k) = c0 * g0(i,j,k) + c1 * g1(i,j,k);
			}			
		}
	}
}

template<typename T, typename R>
class Grid3_XNeumann {
	Grid3<T,R> &grid;

public : 
	typedef typename Grid3<T,R>::index_type index_type;

	Grid3_XNeumann(Grid3<T,R> &_grid) : grid(_grid) {}
	void operator() ( const basis::blocked_range2d<index_type>& r ) const {	
#ifdef USE_OPENMP
		int my_begin = int(r.rows().begin());
		int my_end = int(r.rows().end());
		int i;

#pragma omp parallel for
#else
		index_type my_begin = r.rows().begin();
		index_type my_end = r.rows().end();
		index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			for( index_type j=r.cols().begin(); j!=r.cols().end( ); ++j ) {
				grid(0,i,j) = grid(1,i,j);
				grid(grid.nX-1,i,j) = grid(grid.nX-2,i,j);
			}
		}				
	}	
};

template<typename T, typename R>
class Grid3_YNeumann {
	Grid3<T,R> &grid;

public : 
	typedef typename Grid3<T,R>::index_type index_type;

	Grid3_YNeumann(Grid3<T,R> &_grid) : grid(_grid) {}
	void operator() ( const basis::blocked_range2d<index_type>& r ) const {	
#ifdef USE_OPENMP
		int my_begin = int(r.rows().begin());
		int my_end = int(r.rows().end());
		int i;

#pragma omp parallel for
#else
		index_type my_begin = r.rows().begin();
		index_type my_end = r.rows().end();
		index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			for( index_type j=r.cols().begin(); j!=r.cols().end( ); ++j ) {
				grid(i,0,j) = grid(i,1,j);
				grid(i,grid.nY-1,j) = grid(i,grid.nY-2,j);
			}			
		}										
	}	
};

template<typename T, typename R>
class Grid3_ZNeumann {
	Grid3<T,R> &grid;

public : 
	typedef typename Grid3<T,R>::index_type index_type;

	Grid3_ZNeumann(Grid3<T,R> &_grid) : grid(_grid) {}
	void operator() ( const basis::blocked_range2d<index_type>& r ) const {	
#ifdef USE_OPENMP
		int my_begin = int(r.rows().begin());
		int my_end = int(r.rows().end());
		int i;

#pragma omp parallel for
#else
		index_type my_begin = r.rows().begin();
		index_type my_end = r.rows().end();
		index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			for( index_type j=r.cols().begin(); j!=r.cols().end( ); ++j ) {
				grid(i,j,0) = grid(i,j,1);
				grid(i,j,grid.nZ-1) = grid(i,j,grid.nZ-2);
			}			
		}										
	}	
};

template<typename T, typename R>
class Grid3_XNeumann_Scaling {
	Grid3<T,R> &grid;
	T factor;

public : 
	typedef typename Grid3<T,R>::index_type index_type;

	Grid3_XNeumann_Scaling(Grid3<T,R> &_grid, T _fac) : grid(_grid), factor(_fac) {}
	void operator() ( const basis::blocked_range2d<index_type>& r ) const {	
#ifdef USE_OPENMP
		int my_begin = int(r.rows().begin());
		int my_end = int(r.rows().end());
		int i;

#pragma omp parallel for
#else
		index_type my_begin = r.rows().begin();
		index_type my_end = r.rows().end();
		index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			for( index_type j=r.cols().begin(); j!=r.cols().end( ); ++j ) {
				grid(0,i,j) = grid(1,i,j) * factor;
				grid(grid.nX-1,i,j) = grid(grid.nX-2,i,j) * factor;
			}
		}	
	}	
};

template<typename T, typename R>
class Grid3_YNeumann_Scaling {
	Grid3<T,R> &grid;
	T factor;

public : 
	typedef typename Grid3<T,R>::index_type index_type;

	Grid3_YNeumann_Scaling(Grid3<T,R> &_grid, T _fac) : grid(_grid), factor(_fac) {}
	void operator() ( const basis::blocked_range2d<index_type>& r ) const {	
#ifdef USE_OPENMP
		int my_begin = int(r.rows().begin());
		int my_end = int(r.rows().end());
		int i;

#pragma omp parallel for
#else
		index_type my_begin = r.rows().begin();
		index_type my_end = r.rows().end();
		index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			for( index_type j=r.cols().begin(); j!=r.cols().end( ); ++j ) {
				grid(i,0,j) = grid(i,1,j) * factor;
				grid(i,grid.nY-1,j) = grid(i,grid.nY-2,j) * factor;
			}
		}
	}	
};

template<typename T, typename R>
class Grid3_ZNeumann_Scaling {
	Grid3<T,R> &grid;
	T factor;

public : 
	typedef typename Grid3<T,R>::index_type index_type;

	Grid3_ZNeumann_Scaling(Grid3<T,R> &_grid, T _fac) : grid(_grid), factor(_fac) {}
	void operator() ( const basis::blocked_range2d<index_type>& r ) const {	
#ifdef USE_OPENMP
		int my_begin = int(r.rows().begin());
		int my_end = int(r.rows().end());
		int i;

#pragma omp parallel for
#else
		index_type my_begin = r.rows().begin();
		index_type my_end = r.rows().end();
		index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			for( index_type j=r.cols().begin(); j!=r.cols().end( ); ++j ) {
				grid(i,j,0) = grid(i,j,1) * factor;
				grid(i,j,grid.nZ-1) = grid(i,j,grid.nZ-2) * factor;
			}			
		}										
	}	
};

template<typename T, typename R>
class Grid3_XDirichlet {
	Grid3<T,R> &grid;
	T val;

public : 
	typedef typename Grid3<T,R>::index_type index_type;

	Grid3_XDirichlet(Grid3<T,R> &_grid, T _val) : grid(_grid), val(_val) {}
	void operator() ( const basis::blocked_range2d<index_type>& r ) const {	
#ifdef USE_OPENMP
		int my_begin = int(r.rows().begin());
		int my_end = int(r.rows().end());
		int i;

#pragma omp parallel for
#else
		index_type my_begin = r.rows().begin();
		index_type my_end = r.rows().end();
		index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			for( index_type j=r.cols().begin(); j!=r.cols().end( ); ++j ) {
				grid(0,i,j) = grid(grid.nX-1,i,j) = val;
			}
		}
	}	
};

template<typename T, typename R>
class Grid3_YDirichlet {
	Grid3<T,R> &grid;
	T val;

public : 
	typedef typename Grid3<T,R>::index_type index_type;

	Grid3_YDirichlet(Grid3<T,R> &_grid, T _val) : grid(_grid), val(_val) {}
	void operator() ( const basis::blocked_range2d<index_type>& r ) const {	
#ifdef USE_OPENMP
		int my_begin = int(r.rows().begin());
		int my_end = int(r.rows().end());
		int i;

#pragma omp parallel for
#else
		index_type my_begin = r.rows().begin();
		index_type my_end = r.rows().end();
		index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			for( index_type j=r.cols().begin(); j!=r.cols().end( ); ++j ) {
				grid(i,0,j) = grid(i,grid.nY-1,j) = val;
			}
		}
	}	
};

template<typename T, typename R>
class Grid3_ZDirichlet {
	Grid3<T,R> &grid;
	T val;

public : 
	typedef typename Grid3<T,R>::index_type index_type;

	Grid3_ZDirichlet(Grid3<T,R> &_grid, T _val) : grid(_grid), val(_val) {}
	void operator() ( const basis::blocked_range2d<index_type>& r ) const {	
#ifdef USE_OPENMP
		int my_begin = int(r.rows().begin());
		int my_end = int(r.rows().end());
		int i;

#pragma omp parallel for
#else
		index_type my_begin = r.rows().begin();
		index_type my_end = r.rows().end();
		index_type i;
#endif
		for( i = my_begin; i < my_end; ++i ) {
			for( index_type j=r.cols().begin(); j!=r.cols().end( ); ++j ) {
				grid(i,j,0) = grid(i,j,grid.nZ-1) = val;
			}			
		}										
	}	
};

};		// arithmetic