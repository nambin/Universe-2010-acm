#pragma once

#include <cmath>
#include <vector>
#include "../MathLibrary/MultiGrid.h"
#include "Grid3_PoissonSolver.h"

namespace arithmetic {

template<typename T>
class Grid3_PoissonSolver_GMG : public Grid3_PoissonSolver<T>, public MultiGrid {	
public : 
	friend class MultiGrid;

	enum GRID_TYPE { NODE, CELL };

	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Grid3_PoissonSolver_GMG();
	Grid3_PoissonSolver_GMG(index_type nX, index_type nY, index_type nZ, T h, const Vector3<T>& offset, size_t num_level = 10, GRID_TYPE type = CELL);
	~Grid3_PoissonSolver_GMG() { release(); }

	void release();

	void make_coarse_poisson_solvers(size_t num_level, GRID_TYPE type);

	void set_iteration_number_fine_grid(size_t _iter) { iter_fine_grid = _iter; }
	void set_iteration_number(size_t pre, size_t post) { pre_smoothing = pre; post_smoothing = post; iter_fine_grid = MAX(pre_smoothing, iter_fine_grid); }
	void set_parallel_threshold(index_type thres) { parallel_thres = thres; }

	void scale(const T s);

	void printResidualInfo(index_type, index_type, index_type) const;
	
	//////////////////////////////////////////////////////////////////////////
	// Solver
	const size_t solve(size_t maxItr, T tolerance);	

	static const std::string XmlElementName() { return std::string("Grid3_PoissonSolver_GMG"); }

protected : 
	const size_t max_level() const { return num_level - 1; }
	void set_zero(const size_t lv) { solvers[lv]->x.zero(); }

	void smooth_ith_level(size_t lv, const size_t iter);		// red-black gauss-seidal for level lv

	void restrict_residual(size_t lv);							// restrict residual from level (lv) to level (lv+1)
	void prolongate_correction(size_t lv);						// prolongate solution from level (lv) to level (lv-1)

	const T get_residual(const size_t lv) const;
	void get_residual(Grid3<T,T>& res, const Grid3_PoissonMatrix<T>&, const Grid3<T,T>& x, const Grid3<T,T>& b) const;

	const T get_max_residual(size_t i);
	const T get_l2_residual(size_t i);

public : 
	index_type								parallel_thres;

	size_t									iter_fine_grid;
	size_t									pre_smoothing, post_smoothing;	

	GRID_TYPE								grid_type;		// Now, works for only CELL
	size_t									num_level;

	std::vector<Grid3_PoissonSolver<T>*>	solvers;
	mutable std::vector<Grid3<T,T>>			tmp;
};

template<typename T>
Grid3_PoissonSolver_GMG<T>::Grid3_PoissonSolver_GMG(index_type nX, index_type nY, index_type nZ, T h, const Vector3<T>& offset, size_t num_level, GRID_TYPE type) 
: Grid3_PoissonSolver<T>(nX, nY, nZ, h, offset), iter_fine_grid(2), pre_smoothing(2), post_smoothing(2), parallel_thres(30) {
	make_coarse_poisson_solvers(num_level, type);
}

template<typename T>
void Grid3_PoissonSolver_GMG<T>::release() {	
	for(size_t i=1;i<solvers.size();++i)
		delete solvers[i];
	solvers.clear();
	tmp.clear();

	num_level = 0;
}

template<typename T>
void Grid3_PoissonSolver_GMG<T>::make_coarse_poisson_solvers(size_t num_lv, GRID_TYPE type) {

	release();	

	size_t max_level = 1;
	index_type nX_2 = x.nX-2, nY_2 = x.nY-2, nZ_2 = x.nZ-2;
	const T h = x.dx;

	while(true) {
		if(nX_2 % 2 == 0 && nY_2 % 2 == 0 && nZ_2 % 2 == 0 && nX_2 > 3 && nY_2 > 3 && nZ_2 > 3) {
			max_level++;

			nX_2 /= 2;
			nY_2 /= 2;
			nZ_2 /= 2;
		}
		else 
			break;
	}

	num_level = MIN(max_level, num_lv);		
	grid_type = type;

	solvers.resize(num_level);
	tmp.resize(num_level);

	solvers[0] = this;
	tmp[0].resize(solvers[0]->x);

	int grid_divisor = 1;
	nX_2 = x.nX-2, nY_2 = x.nY-2, nZ_2 = x.nZ-2;

	for(size_t i=1;i<solvers.size();++i) {
		grid_divisor *= 2;
		solvers[i] = new Grid3_PoissonSolver(nX_2 / grid_divisor + 2, nY_2 / grid_divisor + 2, nZ_2 / grid_divisor + 2, h * grid_divisor, x.offset);
		tmp[i].resize(solvers[i]->x);
	}

	log_info << "Grid3_PoissonSolver_GMG _ Multigrid Level _ " << num_level << endl;	
	log_info << "Resolution of Coarsest Grid _ " << nX_2 / grid_divisor << "x" << nY_2 / grid_divisor << "x" << nZ_2 / grid_divisor << endl;
}

template<typename T>
const size_t Grid3_PoissonSolver_GMG<T>::solve(size_t maxItr, T tolerance) {
	const size_t iter = MultiGrid::v_cycle_fully_converging(*this, 0, tolerance, maxItr);		
	return iter;
}

template<typename T>
const T Grid3_PoissonSolver_GMG<T>::get_residual(const size_t lv) const {
	get_residual(tmp[lv], solvers[lv]->poissonMatrix, solvers[lv]->x, solvers[lv]->b);
	return tmp[lv].infNorm();
}

template<typename T>
void Grid3_PoissonSolver_GMG<T>::get_residual(Grid3<T,T>& res, const Grid3_PoissonMatrix<T>& poissonMatrix, const Grid3<T,T>& x, const Grid3<T,T>& b) const {

	const Grid3x4<T>&		poisson = poissonMatrix.poisson;
	const Grid3<bool,T>&	isBoundary = poissonMatrix.isBoundary;

	res.zero();

#ifdef USE_OPENMP
	if(x.nX > parallel_thres) {
		int my_end = int(x.nX-1);
		int i;

#pragma omp parallel for
		for(i=1;i<my_end;++i) {	
			for(index_type j=1;j<x.nY-1;++j) {
				for(index_type k=1;k<x.nZ-1;++k) {
					if(!isBoundary(i,j,k)) {
						res(i,j,k) =	( b(i,j,k)
										- poisson(i-1,j,k,1) * x(i-1,j,k)
										- poisson(i,j,k,1)   * x(i+1,j,k)
										- poisson(i,j-1,k,2) * x(i,j-1,k)
										- poisson(i,j,k,2)   * x(i,j+1,k) 
										- poisson(i,j,k-1,3) * x(i,j,k-1)
										- poisson(i,j,k,3)   * x(i,j,k+1) 

										- poisson(i,j,k,0)   * x(i,j,k) ) ;	
					}
				}
			}
		}
	}
	else 
#endif

	{
		index_type my_end = x.nX-1, i;

		for(i=1;i<my_end;++i) {	
			for(index_type j=1;j<x.nY-1;++j) {
				for(index_type k=1;k<x.nZ-1;++k) {
					if(!isBoundary(i,j,k)) {
						res(i,j,k) =	( b(i,j,k)
										- poisson(i-1,j,k,1) * x(i-1,j,k)
										- poisson(i,j,k,1)   * x(i+1,j,k)
										- poisson(i,j-1,k,2) * x(i,j-1,k)
										- poisson(i,j,k,2)   * x(i,j+1,k) 
										- poisson(i,j,k-1,3) * x(i,j,k-1)
										- poisson(i,j,k,3)   * x(i,j,k+1) 

										- poisson(i,j,k,0)   * x(i,j,k) ) ;	
					}
				}
			}
		}
	}
}

template<typename T>
void Grid3_PoissonSolver_GMG<T>::printResidualInfo(index_type i, index_type j, index_type k) const {
	log_info << i << " " << j << " " << k << " Residual info " << endl;

	const Grid3x4<T>& poisson = poissonMatrix.poisson;

	log_info << "b : " << b(i,j,k) << endl;
	log_info << "x -1 _ : " << poisson(i-1,j,k,1)	<< " " << x(i-1,j,k) << " " << + (x(i,j,k) - x(i-1,j,k)) * b.dxInv << endl;
	log_info << "x +1 _ : " << poisson(i  ,j,k,1)   << " " << x(i+1,j,k) << " " << - (x(i,j,k) - x(i+1,j,k)) * b.dxInv << endl;
	log_info << "y -1 _ : " << poisson(i,j-1,k,2)	<< " " << x(i,j-1,k) << " " << + (x(i,j,k) - x(i,j-1,k)) * b.dxInv << endl;
	log_info << "y +1 _ : " << poisson(i,j  ,k,2)   << " " << x(i,j+1,k) << " " << - (x(i,j,k) - x(i,j+1,k)) * b.dxInv << endl;
	log_info << "z -1 _ : " << poisson(i,j,k-1,3)	<< " " << x(i,j,k-1) << " " << + (x(i,j,k) - x(i,j,k-1)) * b.dxInv << endl;
	log_info << "z +1 _ : " << poisson(i,j,k  ,3)   << " " << x(i,j,k+1) << " " << - (x(i,j,k) - x(i,j,k+1)) * b.dxInv << endl;
	log_info << "center : " << poisson(i,j,k,0)     << " " << x(i,j,k) << endl;
	

	T res = b(i,j,k)
			- poisson(i-1,j,k,1) * x(i-1,j,k)
			- poisson(i,j,k,1)   * x(i+1,j,k)
			- poisson(i,j-1,k,2) * x(i,j-1,k)
			- poisson(i,j,k,2)   * x(i,j+1,k) 
			- poisson(i,j,k-1,3) * x(i,j,k-1)
			- poisson(i,j,k,3)   * x(i,j,k+1) 

			- poisson(i,j,k,0)   * x(i,j,k) ;	

	log_info << "Residual : " << res << endl;
}

template<typename T>
void Grid3_PoissonSolver_GMG<T>::smooth_ith_level(size_t ith, const size_t _iter) {

	const size_t num_iter = _iter;

	const Grid3x4<T>&		poisson = solvers[ith]->poissonMatrix.poisson;
	const Grid3<bool,T>&	isBoundary = solvers[ith]->poissonMatrix.isBoundary;

	const Grid3<T,T>&		b = solvers[ith]->b;
	Grid3<T,T>&				x = solvers[ith]->x;
	Grid3<T,T>&				x_tmp = tmp[ith];

	x_tmp.zero();
	for(size_t iter = 0; iter < num_iter ; ++iter) {

#ifdef USE_OPENMP
		if(x.nX > parallel_thres) {
			int my_end = int(x.nX-1), i;

#pragma omp parallel for
			for( i = 1; i < my_end; ++i ) {		
				for( index_type j = 1; j < x.nY-1; ++j) {

					index_type start_k;

					if( (i+j) % 2 == 0) start_k = 1;
					else				start_k = 2;

					for(index_type k=start_k;k<x.nZ-1;k+=2) {
						if(!isBoundary(i,j,k)) {
							x_tmp(i,j,k) =	( b(i,j,k)
											- poisson(i-1,j,k,1) * x(i-1,j,k)
											- poisson(i,j,k,1)   * x(i+1,j,k)
											- poisson(i,j-1,k,2) * x(i,j-1,k)
											- poisson(i,j,k,2)   * x(i,j+1,k) 
											- poisson(i,j,k-1,3) * x(i,j,k-1)
											- poisson(i,j,k,3)   * x(i,j,k+1) ) / poisson(i,j,k,0);			
						}
					}

				}				
			}

			my_end = int(x.nX-1);

#pragma omp parallel for
			for( i = 1; i < my_end; ++i ) {
				for( index_type j = 1; j < x.nY-1; ++j) {

					index_type start_k;

					if( (i+j) % 2 == 0) start_k = 2;
					else				start_k = 1;

					for(index_type k=start_k;k<x.nZ-1;k+=2) {
						if(!isBoundary(i,j,k)) {
							x_tmp(i,j,k) =	( b(i,j,k)
											- poisson(i-1,j,k,1) * x_tmp(i-1,j,k)
											- poisson(i,j,k,1)   * x_tmp(i+1,j,k)
											- poisson(i,j-1,k,2) * x_tmp(i,j-1,k)
											- poisson(i,j,k,2)   * x_tmp(i,j+1,k) 
											- poisson(i,j,k-1,3) * x_tmp(i,j,k-1)
											- poisson(i,j,k,3)   * x_tmp(i,j,k+1) ) / poisson(i,j,k,0);			
						}
					}

				}
			}
		}
		else 
#endif
		{	
			index_type my_end = x.nX-1, i;			
			for( i = 1; i < my_end; ++i ) {		
				for( index_type j = 1; j < x.nY-1; ++j) {

					index_type start_k;

					if( (i+j) % 2 == 0) start_k = 1;
					else				start_k = 2;

					for(index_type k=start_k;k<x.nZ-1;k+=2) {
						if(!isBoundary(i,j,k)) {
							x_tmp(i,j,k) =	( b(i,j,k)
											- poisson(i-1,j,k,1) * x(i-1,j,k)
											- poisson(i,j,k,1)   * x(i+1,j,k)
											- poisson(i,j-1,k,2) * x(i,j-1,k)
											- poisson(i,j,k,2)   * x(i,j+1,k) 
											- poisson(i,j,k-1,3) * x(i,j,k-1)
											- poisson(i,j,k,3)   * x(i,j,k+1) ) / poisson(i,j,k,0);			
						}
					}

				}	
			}

			my_end = x.nX-1;
			for( i = 1; i < my_end; ++i ) {
				for( index_type j = 1; j < x.nY-1; ++j) {

					index_type start_k;

					if( (i+j) % 2 == 0) start_k = 2;
					else				start_k = 1;

					for(index_type k=start_k;k<x.nZ-1;k+=2) {
						if(!isBoundary(i,j,k)) {
							x_tmp(i,j,k) =	( b(i,j,k)
											- poisson(i-1,j,k,1) * x_tmp(i-1,j,k)
											- poisson(i,j,k,1)   * x_tmp(i+1,j,k)
											- poisson(i,j-1,k,2) * x_tmp(i,j-1,k)
											- poisson(i,j,k,2)   * x_tmp(i,j+1,k) 
											- poisson(i,j,k-1,3) * x_tmp(i,j,k-1)
											- poisson(i,j,k,3)   * x_tmp(i,j,k+1) ) / poisson(i,j,k,0);			
						}
					}
				}
			}

			//index_type my_end = x.nX-1, i;			
			//for( i = 1; i < my_end; ++i ) {		
			//	for( index_type j = 1; j < x.nY-1; ++j) {
			//		for(index_type k=1;k<x.nZ-1; ++k) {
			//			if(!isBoundary(i,j,k)) {
			//				x_tmp(i,j,k) =	( b(i,j,k)
			//								- poisson(i-1,j,k,1) * x_tmp(i-1,j,k)
			//								- poisson(i,j,k,1)   * x(i+1,j,k)
			//								- poisson(i,j-1,k,2) * x_tmp(i,j-1,k)
			//								- poisson(i,j,k,2)   * x(i,j+1,k) 
			//								- poisson(i,j,k-1,3) * x_tmp(i,j,k-1)
			//								- poisson(i,j,k,3)   * x(i,j,k+1) ) / poisson(i,j,k,0);			
			//			}
			//		}

			//	}	
			//}

		}

		swap(x_tmp, x);
	}
}

template<typename T>
void Grid3_PoissonSolver_GMG<T>::restrict_residual(size_t ith) {

	const real eighth = real(1.0/8.0);

	Grid3<T,T>&	res = tmp[ith];
	Grid3<T,T>&	res_coarse = solvers[ith+1]->b;

	get_residual(res, solvers[ith]->poissonMatrix, solvers[ith]->x, solvers[ith]->b);

	res_coarse.zero();

#ifdef USE_OPENMP		
	if(res_coarse.nX > parallel_thres) {
		int my_end = int(res_coarse.nX-1), i;

#pragma omp parallel for
		for(i=1;i<my_end;++i) 
			for(index_type j=1;j<res_coarse.nY-1;++j) 
				for(index_type k=1;k<res_coarse.nZ-1;++k) 
					res_coarse(i,j,k) = (res(i*2-1, j*2-1, k*2-1) + res(i*2-1, j*2, k*2-1) + res(i*2, j*2-1, k*2-1) + res(i*2, j*2, k*2-1)
									   + res(i*2-1, j*2-1, k*2  ) + res(i*2-1, j*2, k*2  ) + res(i*2, j*2-1, k*2  ) + res(i*2, j*2, k*2  )) * eighth;

	}
	else 
#endif
	{
		index_type my_end = res_coarse.nX-1, i;
		for(i=1;i<my_end;++i) 
			for(index_type j=1;j<res_coarse.nY-1;++j) 
				for(index_type k=1;k<res_coarse.nZ-1;++k) 
					res_coarse(i,j,k) = (res(i*2-1, j*2-1, k*2-1) + res(i*2-1, j*2, k*2-1) + res(i*2, j*2-1, k*2-1) + res(i*2, j*2, k*2-1)
									   + res(i*2-1, j*2-1, k*2  ) + res(i*2-1, j*2, k*2  ) + res(i*2, j*2-1, k*2  ) + res(i*2, j*2, k*2  ) ) * eighth;
	}
}

template<typename T>
void Grid3_PoissonSolver_GMG<T>::prolongate_correction(size_t ith) {	
	const Grid3<bool,T>&	isbdry = solvers[ith-1]->poissonMatrix.isBoundary;
	const Grid3<T,T>&	x = solvers[ith]->x;	
	Grid3<T,T>&			prolong_x = solvers[ith-1]->x;

#ifdef USE_OPENMP		
	if(x.nX > parallel_thres) {
		int my_end = int(x.nX-1), i;

#pragma omp parallel for
		for(i=1;i<my_end;++i) {	
			for(index_type j=1;j<x.nY-1;++j) {
				for(index_type k=1;k<x.nZ-1;++k) {
					if(!isbdry(i*2-1,j*2-1,k*2-1)) prolong_x(i*2-1, j*2-1, k*2-1) += (x(i,j,k));
					if(!isbdry(i*2-1,j*2  ,k*2-1)) prolong_x(i*2-1, j*2  , k*2-1) += (x(i,j,k));
					if(!isbdry(i*2  ,j*2-1,k*2-1)) prolong_x(i*2  , j*2-1, k*2-1) += (x(i,j,k));
					if(!isbdry(i*2  ,j*2  ,k*2-1)) prolong_x(i*2  , j*2  , k*2-1) += (x(i,j,k));
					if(!isbdry(i*2-1,j*2-1,k*2  )) prolong_x(i*2-1, j*2-1, k*2  ) += (x(i,j,k));
					if(!isbdry(i*2-1,j*2  ,k*2  )) prolong_x(i*2-1, j*2  , k*2  ) += (x(i,j,k));
					if(!isbdry(i*2  ,j*2-1,k*2  )) prolong_x(i*2  , j*2-1, k*2  ) += (x(i,j,k));
					if(!isbdry(i*2  ,j*2  ,k*2  )) prolong_x(i*2  , j*2  , k*2  ) += (x(i,j,k));
				}
			}
		}

	}
	else 
#endif
	{
		index_type my_end = x.nX-1, i;
		for(i=1;i<my_end;++i) {	
			for(index_type j=1;j<x.nY-1;++j) {
				for(index_type k=1;k<x.nZ-1;++k) {
					if(!isbdry(i*2-1,j*2-1,k*2-1)) prolong_x(i*2-1, j*2-1, k*2-1) += (x(i,j,k));
					if(!isbdry(i*2-1,j*2  ,k*2-1)) prolong_x(i*2-1, j*2  , k*2-1) += (x(i,j,k));
					if(!isbdry(i*2  ,j*2-1,k*2-1)) prolong_x(i*2  , j*2-1, k*2-1) += (x(i,j,k));
					if(!isbdry(i*2  ,j*2  ,k*2-1)) prolong_x(i*2  , j*2  , k*2-1) += (x(i,j,k));
					if(!isbdry(i*2-1,j*2-1,k*2  )) prolong_x(i*2-1, j*2-1, k*2  ) += (x(i,j,k));
					if(!isbdry(i*2-1,j*2  ,k*2  )) prolong_x(i*2-1, j*2  , k*2  ) += (x(i,j,k));
					if(!isbdry(i*2  ,j*2-1,k*2  )) prolong_x(i*2  , j*2-1, k*2  ) += (x(i,j,k));
					if(!isbdry(i*2  ,j*2  ,k*2  )) prolong_x(i*2  , j*2  , k*2  ) += (x(i,j,k));
				}
			}
		}
	}
}

template<typename T>
const T Grid3_PoissonSolver_GMG<T>::get_max_residual(size_t ith) {
	const Grid3x4<T>&		poisson = solvers[ith]->poissonMatrix.poisson;
	const Grid3<bool,T>&	isBoundary = solvers[ith]->poissonMatrix.isBoundary;
	const Grid3<T,T>&		x = solvers[ith]->x;
	const Grid3<T,T>&		b = solvers[ith]->b;

	T max_res = 0;	index_type i_max, j_max, k_max;
	for(index_type i=1;i<x.nX-1;++i) {
		for(index_type j=1;j<x.nY-1;++j) {
			for(index_type k=1;k<x.nZ-1;++k) {
				if(!isBoundary(i,j,k)) {
					real res =	( b(i,j,k)
								- poisson(i-1,j,k,1) * x(i-1,j,k)
								- poisson(i,j,k,1)   * x(i+1,j,k)
								- poisson(i,j-1,k,2) * x(i,j-1,k)
								- poisson(i,j,k,2)   * x(i,j+1,k) 
								- poisson(i,j,k-1,3) * x(i,j,k-1)
								- poisson(i,j,k,3)   * x(i,j,k+1) 

								- poisson(i,j,k,0)   * x(i,j,k) ) ;	

					if(max_res < ::abs(res)) {
						max_res = ::abs(res);
						i_max = i; j_max = j; k_max = k;
					}
				}				
			}			
		}
	}

	log_info << "Max Residual for level " << ith << " _ Grid3_PoissonSolver_GMG : " << max_res << " " << i_max << " " << j_max << " " << k_max << endl;

	return max_res;
}

template<typename T>
const T Grid3_PoissonSolver_GMG<T>::get_l2_residual(size_t ith) {
	const Grid3x4<T>&		poisson = solvers[ith]->poissonMatrix.poisson;
	const Grid3<bool,T>&	isBoundary = solvers[ith]->poissonMatrix.isBoundary;
	const Grid3<T,T>&		x = solvers[ith]->x;
	const Grid3<T,T>&		b = solvers[ith]->b;

	T norml2 = 0; size_t n_s = 0;
	for(index_type i=1;i<x.nX-1;++i) {
		for(index_type j=1;j<x.nY-1;++j) {
			for(index_type k=1;k<x.nZ-1;++k) {
				if(!isBoundary(i,j,k)) {
					real res =	( b(i,j,k)
								- poisson(i-1,j,k,1) * x(i-1,j,k)
								- poisson(i,j,k,1)   * x(i+1,j,k)
								- poisson(i,j-1,k,2) * x(i,j-1,k)
								- poisson(i,j,k,2)   * x(i,j+1,k) 
								- poisson(i,j,k-1,3) * x(i,j,k-1)
								- poisson(i,j,k,3)   * x(i,j,k+1) 

								- poisson(i,j,k,0)   * x(i,j,k) ) ;	

					norml2 += (res * res);
					n_s++;
				}				
			}	
		}
	}

	log_info << "L2 Residual for level " << ith << " _ Grid3_PoissonSolver_GMG : " << sqrt(norml2 / n_s) << endl;

	return sqrt(norml2 / n_s);
}

template<typename T>
void Grid3_PoissonSolver_GMG<T>::scale(const T s) {
	poissonMatrix.poisson *= s;
	b *= s;

	for(size_t i=1;i<solvers.size();++i)
		solvers[i]->scale(s);
}

};		// arithmetic