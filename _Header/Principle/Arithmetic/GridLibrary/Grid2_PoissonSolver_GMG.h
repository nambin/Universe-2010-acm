#pragma once

#include <cmath>
#include <vector>
#include "../MathLibrary/MultiGrid.h"
#include "Grid2_PoissonSolver.h"

namespace arithmetic {

template<typename T>
class Grid2_PoissonSolver_GMG : public Grid2_PoissonSolver<T>, public MultiGrid {	
public : 
	friend class MultiGrid; 

	enum GRID_TYPE { NODE, CELL };

	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Grid2_PoissonSolver_GMG();
	Grid2_PoissonSolver_GMG(index_type nX, index_type nY, T h, const Vector2<T>& offset, size_t num_level = 15, GRID_TYPE type = CELL);
	~Grid2_PoissonSolver_GMG() { release(); }
	
	void release();

	void make_coarse_poisson_solvers(size_t num_level, GRID_TYPE type);

	void set_iteration_number_fine_grid(size_t _iter) { iter_fine_grid = _iter; }
	void set_iteration_number(size_t pre, size_t post) { pre_smoothing = pre; post_smoothing = post; iter_fine_grid = MAX(pre_smoothing, iter_fine_grid); }
	void set_parallel_threshold(index_type thres) { parallel_thres = thres; }

	void scale(const T s);

	//////////////////////////////////////////////////////////////////////////
	// Solver
	const size_t solve(size_t maxItr, T tolerance);	

	const bool importBinary(std::ifstream&);
	const bool exportBinary(std::ofstream&) const;

	static const std::string XmlElementName() { return std::string("Grid2_PoissonSolver_GMG"); }

protected : 
	const size_t max_level() const { return num_level - 1; }
	void set_zero(const size_t lv) { solvers[lv]->x.zero(); }

	void smooth_ith_level(size_t lv, const size_t iter);			// red-black gauss-seidal for level lv
	
	void restrict_residual(size_t lv);								// restrict residual from level (lv) to level (lv+1)
	void prolongate_correction(size_t lv);							// prolongate solution from level (lv) to level (lv-1)

	const T get_residual(const size_t lv) const;
	void get_residual(Grid2<T,T>& res, const Grid2_PoissonMatrix<T>&, const Grid2<T,T>& x, const Grid2<T,T>& b) const;

	const T get_max_residual(size_t i);
	const T get_l2_residual(size_t i);

public : 
	index_type								parallel_thres;

	size_t									iter_fine_grid;
	size_t									pre_smoothing, post_smoothing;	

	GRID_TYPE								grid_type;		// Now, works for only CELL
	size_t									num_level;

	std::vector<Grid2_PoissonSolver<T>*>	solvers;
	mutable std::vector<Grid2<T,T>>			tmp;
};

template<typename T>
Grid2_PoissonSolver_GMG<T>::Grid2_PoissonSolver_GMG()
: iter_fine_grid(2), pre_smoothing(2), post_smoothing(2), parallel_thres(30) 
{}

template<typename T>
Grid2_PoissonSolver_GMG<T>::Grid2_PoissonSolver_GMG(index_type nX, index_type nY, T h, const Vector2<T>& offset, size_t num_level, GRID_TYPE type) 
: Grid2_PoissonSolver<T>(nX, nY, h, offset), iter_fine_grid(2), pre_smoothing(2), post_smoothing(2), parallel_thres(30) {
	make_coarse_poisson_solvers(num_level, type);
}

template<typename T>
void Grid2_PoissonSolver_GMG<T>::release() {	
	for(size_t i=1;i<solvers.size();++i)
		delete solvers[i];
	solvers.clear();
	tmp.clear();

	num_level = 0;
}

template<typename T>
void Grid2_PoissonSolver_GMG<T>::make_coarse_poisson_solvers(size_t num_lv, GRID_TYPE type) {

	release();	

	size_t max_level = 1;
	index_type nX_2 = b.nX-2, nY_2 = b.nY-2;
	const T h = b.dx;

	while(true) {
		if(nX_2 % 2 == 0 && nY_2 % 2 == 0 && nX_2 > 3 && nY_2 > 3) {
			max_level++;

			nX_2 /= 2;
			nY_2 /= 2;
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
	nX_2 = b.nX-2, nY_2 = b.nY-2;

	for(size_t i=1;i<solvers.size();++i) {
		grid_divisor *= 2;
		solvers[i] = new Grid2_PoissonSolver(nX_2 / grid_divisor + 2, nY_2 / grid_divisor + 2, h * grid_divisor, x.offset);
		tmp[i].resize(solvers[i]->x);
	}

	log_info << "Grid2_PoissonSolver_GMG _ Multigrid Level _ " << num_level << endl;
	log_info << "Resolution of Coarsest Grid _ " << nX_2 / grid_divisor << "x" << nY_2 / grid_divisor << endl;
}

template<typename T>
const size_t Grid2_PoissonSolver_GMG<T>::solve(size_t maxItr, T tolerance) {

	//restrict_residual(0);
	//smooth_ith_level(0, 10000);		get_max_residual(0);
	//smooth_ith_level(1, 3000);		get_max_residual(1);

	//T sum = 0;
	//for(index_type i=1;i<solvers[1]->x.nX-1;++i) {
	//	for(index_type j=1;j<solvers[1]->x.nY-1;++j) {
	//		const T avg = (solvers[0]->x(i*2-1,j*2-1) + solvers[0]->x(i*2-1,j*2) + solvers[0]->x(i*2,j*2-1) + solvers[0]->x(i*2,j*2)) * 0.25f;
	//		sum += ::abs(avg - solvers[1]->x(i,j));
	//	}
	//}
	//log_info << "L1 Error : " << sum / ((solvers[1]->x.nX-2) * (solvers[1]->x.nY-2)) << endl;

	//index_type i,j;

	//i = 1; j = 1;
	//log_info << solvers[1]->x(i,j) << endl;
	//log_info << solvers[0]->x(i*2-1,j*2-1) << " " << solvers[0]->x(i*2-1,j*2) << " " << solvers[0]->x(i*2,j*2-1) << " " << solvers[0]->x(i*2,j*2) << endl;

	//i = 2; j = 4;
	//log_info << solvers[1]->x(i,j) << endl;
	//log_info << solvers[0]->x(i*2-1,j*2-1) << " " << solvers[0]->x(i*2-1,j*2) << " " << solvers[0]->x(i*2,j*2-1) << " " << solvers[0]->x(i*2,j*2) << endl;

	//i = 16; j = 10;
	//log_info << solvers[1]->x(i,j) << endl;
	//log_info << solvers[0]->x(i*2-1,j*2-1) << " " << solvers[0]->x(i*2-1,j*2) << " " << solvers[0]->x(i*2,j*2-1) << " " << solvers[0]->x(i*2,j*2) << endl;

//	return MultiGrid::v_cycle(*this, tolerance, maxItr);
	return MultiGrid::v_cycle_fully_converging(*this, 0, tolerance, maxItr);	
}

template<typename T>
const T Grid2_PoissonSolver_GMG<T>::get_residual(const size_t lv) const {
	get_residual(tmp[lv], solvers[lv]->poissonMatrix, solvers[lv]->x, solvers[lv]->b);
	return tmp[lv].infNorm();
}

template<typename T>
void Grid2_PoissonSolver_GMG<T>::get_residual(Grid2<T,T>& res, const Grid2_PoissonMatrix<T>& poissonMatrix, const Grid2<T,T>& x, const Grid2<T,T>& b) const {

	const Grid2x3<T>&		poisson = poissonMatrix.poisson;
	const Grid2<bool,T>&	isBoundary = poissonMatrix.isBoundary;

	res.zero();

#ifdef USE_OPENMP
	if(x.nX > parallel_thres) {
		int my_end = int(x.nX-1);
		int i;

#pragma omp parallel for
		for(i=1;i<my_end;++i) {	
			for(index_type j=1;j<x.nY-1;++j) {
				if(!isBoundary(i,j)) {
					res(i,j) =	( b(i,j)
								- poisson(i-1,j,1) * x(i-1,j)
								- poisson(i,j,1)   * x(i+1,j)
								- poisson(i,j-1,2) * x(i,j-1)
								- poisson(i,j,2)   * x(i,j+1) 

								- poisson(i,j,0)   * x(i,j) ) ;	
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
				if(!isBoundary(i,j)) {
					res(i,j) =	( b(i,j)
								- poisson(i-1,j,1) * x(i-1,j)
								- poisson(i,j,1)   * x(i+1,j)
								- poisson(i,j-1,2) * x(i,j-1)
								- poisson(i,j,2)   * x(i,j+1) 

								- poisson(i,j,0)   * x(i,j) ) ;	
				}
			}
		}
	}
}

template<typename T>
void Grid2_PoissonSolver_GMG<T>::smooth_ith_level(size_t ith, const size_t _iter) {

	const size_t num_iter = _iter;

	const Grid2x3<T>&		poisson = solvers[ith]->poissonMatrix.poisson;
	const Grid2<bool,T>&	isBoundary = solvers[ith]->poissonMatrix.isBoundary;

	const Grid2<T,T>&		b = solvers[ith]->b;
	Grid2<T,T>&				x = solvers[ith]->x;
	Grid2<T,T>&				x_tmp = tmp[ith];

	x_tmp.zero();
	for(size_t iter = 0; iter < num_iter ; ++iter) {

#ifdef USE_OPENMP
		if(x.nX > parallel_thres) {
			int my_end = int(x.nX-1), i;

#pragma omp parallel for
			for( i = 1; i < my_end; ++i ) {		
				index_type start_j;

				if( i % 2 == 0) start_j = 1;
				else			start_j = 2;

				for(index_type j=start_j;j<x.nY-1;j+=2) {
					if(!isBoundary(i,j)) {
						x_tmp(i,j) =	( b(i,j)
										- poisson(i-1,j,1) * x(i-1,j)
										- poisson(i,j,1)   * x(i+1,j)
										- poisson(i,j-1,2) * x(i,j-1)
										- poisson(i,j,2)   * x(i,j+1) ) / poisson(i,j,0);			
					}
				}
			}

			my_end = int(x.nX-1);

#pragma omp parallel for
			for( i = 1; i < my_end; ++i ) {
				index_type start_j;

				if( i % 2 == 0) start_j = 2;
				else			start_j = 1;

				for(index_type j=start_j;j<x.nY-1;j+=2) {
					if(!isBoundary(i,j)) {
						x_tmp(i,j) =	( b(i,j)
										- poisson(i-1,j,1) * x_tmp(i-1,j)
										- poisson(i,j,1)   * x_tmp(i+1,j)
										- poisson(i,j-1,2) * x_tmp(i,j-1)
										- poisson(i,j,2)   * x_tmp(i,j+1) ) / poisson(i,j,0);			
					}
				}
			}
		}
		else 
#endif
		{	
			index_type my_end = x.nX-1, i;			
			for( i = 1; i < my_end; ++i ) {		
				index_type start_j;

				if( i % 2 == 0) start_j = 1;
				else			start_j = 2;

				for(index_type j=start_j;j<x.nY-1;j+=2) {
					if(!isBoundary(i,j)) {
						x_tmp(i,j) =	( b(i,j)
										- poisson(i-1,j,1) * x(i-1,j)
										- poisson(i,j,1)   * x(i+1,j)
										- poisson(i,j-1,2) * x(i,j-1)
										- poisson(i,j,2)   * x(i,j+1) ) / poisson(i,j,0);			
					}
				}
			}

			my_end = x.nX-1;
			for( i = 1; i < my_end; ++i ) {
				index_type start_j;

				if( i % 2 == 0) start_j = 2;
				else			start_j = 1;

				for(index_type j=start_j;j<x.nY-1;j+=2) {
					if(!isBoundary(i,j)) {
						x_tmp(i,j) =	( b(i,j)
										- poisson(i-1,j,1) * x_tmp(i-1,j)
										- poisson(i,j,1)   * x_tmp(i+1,j)
										- poisson(i,j-1,2) * x_tmp(i,j-1)
										- poisson(i,j,2)   * x_tmp(i,j+1) ) / poisson(i,j,0);			
					}
				}
			}
		}

		swap(x_tmp, x);
	}
}

template<typename T>
void Grid2_PoissonSolver_GMG<T>::restrict_residual(size_t ith) {

	const T fourth = T(1.0/4.0);

	Grid2<T,T>&	res = tmp[ith];
	Grid2<T,T>&	res_coarse = solvers[ith+1]->b;

	get_residual(res, solvers[ith]->poissonMatrix, solvers[ith]->x, solvers[ith]->b);

	res_coarse.zero();

#ifdef USE_OPENMP		
	if(res_coarse.nX > parallel_thres) {
		int my_end = int(res_coarse.nX-1), i;

#pragma omp parallel for
		for(i=1;i<my_end;++i) 
			for(index_type j=1;j<res_coarse.nY-1;++j) 
				res_coarse(i,j) = (res(i*2-1, j*2-1) + res(i*2-1, j*2) + res(i*2, j*2-1) + res(i*2, j*2)) * fourth;
	}
	else 
#endif
	{
		index_type my_end = res_coarse.nX-1, i;
		for(i=1;i<my_end;++i) 
			for(index_type j=1;j<res_coarse.nY-1;++j) 
				res_coarse(i,j) = (res(i*2-1, j*2-1) + res(i*2-1, j*2) + res(i*2, j*2-1) + res(i*2, j*2)) * fourth;
	}
}

template<typename T>
void Grid2_PoissonSolver_GMG<T>::prolongate_correction(size_t ith) {		
	const Grid2<bool,T>&	isbdry = solvers[ith-1]->poissonMatrix.isBoundary;
	const Grid2<T,T>&		x = solvers[ith]->x;	
	Grid2<T,T>&				prolong_x = solvers[ith-1]->x;

#ifdef USE_OPENMP		
	if(x.nX > parallel_thres) {
		int my_end = int(x.nX-1), i;

#pragma omp parallel for
		for(i=1;i<my_end;++i) {	
			for(index_type j=1;j<x.nY-1;++j) {

				if(!isbdry(i*2-1, j*2-1))	prolong_x(i*2-1, j*2-1) += (x(i,j));
				if(!isbdry(i*2-1, j*2  ))	prolong_x(i*2-1, j*2  ) += (x(i,j));
				if(!isbdry(i*2  , j*2-1))	prolong_x(i*2  , j*2-1) += (x(i,j));
				if(!isbdry(i*2  , j*2  ))	prolong_x(i*2  , j*2  ) += (x(i,j));

			}
		}
	}
	else 
#endif
	{
		index_type my_end = x.nX-1, i;
		for(i=1;i<my_end;++i) {	
			for(index_type j=1;j<x.nY-1;++j) {

				if(!isbdry(i*2-1, j*2-1))	prolong_x(i*2-1, j*2-1) += (x(i,j));
				if(!isbdry(i*2-1, j*2  ))	prolong_x(i*2-1, j*2  ) += (x(i,j));
				if(!isbdry(i*2  , j*2-1))	prolong_x(i*2  , j*2-1) += (x(i,j));
				if(!isbdry(i*2  , j*2  ))	prolong_x(i*2  , j*2  ) += (x(i,j));

			}
		}
	}
}

template<typename T>
const T Grid2_PoissonSolver_GMG<T>::get_max_residual(size_t ith) {
	const Grid2x3<T>&		poisson = solvers[ith]->poissonMatrix.poisson;
	const Grid2<bool,T>&	isBoundary = solvers[ith]->poissonMatrix.isBoundary;
	const Grid2<T,T>&		x = solvers[ith]->x;
	const Grid2<T,T>&		b = solvers[ith]->b;

	T max_res = 0;	index_type i_max, j_max;
	for(index_type i=1;i<x.nX-1;++i) {
		for(index_type j=1;j<x.nY-1;++j) {
			if(!isBoundary(i,j)) {
				const T res = ( b(i,j)
							- poisson(i-1,j,1) * x(i-1,j)
							- poisson(i,j,1)   * x(i+1,j)
							- poisson(i,j-1,2) * x(i,j-1)
							- poisson(i,j,2)   * x(i,j+1) 

							- poisson(i,j,0)   * x(i,j) ) ;	

				if(max_res < ::abs(res)) {
					max_res = ::abs(res);
					i_max = i; j_max = j;
				}
			}
		}
	}

	log_info << "Max Residual for level " << ith << " _ Grid2_PoissonSolver_GMG : " << max_res << " " << i_max << " " << j_max << endl;

	return max_res;
}

template<typename T>
const T Grid2_PoissonSolver_GMG<T>::get_l2_residual(size_t ith) {
	const Grid2x3<T>&		poisson = solvers[ith]->poissonMatrix.poisson;
	const Grid2<bool,T>&	isBoundary = solvers[ith]->poissonMatrix.isBoundary;
	const Grid2<T,T>&		x = solvers[ith]->x;
	const Grid2<T,T>&		b = solvers[ith]->b;

	T norml2 = 0; size_t n_s = 0;
	for(index_type i=1;i<x.nX-1;++i) {
		for(index_type j=1;j<x.nY-1;++j) {
			if(!isBoundary(i,j)) {
				const T res = ( b(i,j)
							- poisson(i-1,j,1) * x(i-1,j)
							- poisson(i,j,1)   * x(i+1,j)
							- poisson(i,j-1,2) * x(i,j-1)
							- poisson(i,j,2)   * x(i,j+1) 

							- poisson(i,j,0)   * x(i,j) ) ;	

				norml2 += (res * res);
				n_s++;
			}
		}
	}

	log_info << "L2 Residual for level " << ith << " _ Grid2_PoissonSolver_GMG : " << sqrt(norml2 / n_s) << endl;

	return sqrt(norml2 / n_s);
}

template<typename T>
const bool Grid2_PoissonSolver_GMG<T>::importBinary(std::ifstream& ifs) {
	release();

	Grid2_PoissonSolver<T>::importBinary(ifs);

	//	int _type;
	//	ifs.read((char*)(&_type), sizeof(int)); grid_type = GRID_TYPE(_type);

	unsigned int n_lv;
	ifs.read((char*)(&n_lv), sizeof(unsigned int));

	num_level = n_lv;

	solvers.resize(num_level);
	tmp.resize(num_level);

	solvers[0] = this;
	tmp[0].resize(this->x);
	for(size_t i=1;i<solvers.size();++i) {
		solvers[i] = new Grid2_PoissonSolver<T>();
		solvers[i]->importBinary(ifs);

		tmp[i].resize(solvers[i]->x);
	}

	return true;
}

template<typename T>
const bool Grid2_PoissonSolver_GMG<T>::exportBinary(std::ofstream& ofs) const {
	Grid2_PoissonSolver<T>::exportBinary(ofs);

	//	const int _type = int(grid_type);
	//	ofs.write((char*)(&_type), sizeof(int));

	unsigned int n_lv = unsigned int(num_level);
	ofs.write((char*)(&n_lv), sizeof(unsigned int));

	for(size_t i=1;i<solvers.size();++i)
		solvers[i]->exportBinary(ofs);

	return true;
}

template<typename T>
void Grid2_PoissonSolver_GMG<T>::scale(const T s) {
	poissonMatrix.poisson *= s;
	b *= s;

	for(size_t i=1;i<solvers.size();++i)
		solvers[i]->scale(s);
}

};		// arithmetic