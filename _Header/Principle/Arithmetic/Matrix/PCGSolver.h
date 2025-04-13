#ifndef _ARITHMETIC_PCG_SOLVER_H_
#define _ARITHMETIC_PCG_SOLVER_H_

#include <cmath>
#include "FixedSparseMatrix.h"

namespace arithmetic {

template<class T> class SparseColumnLowerFactor {
public : 
	explicit SparseColumnLowerFactor(size_t _n = 0);

	void clear(void);
	void resize(size_t _n);

	void makePreconditioner_IC0(const SparseMatrix<T>&, T modification_parameter = 0.97f, T min_diagonal_ratio = 0.25f);
	void applyPreconditioner(VectorN<T> &ret, const VectorN<T> &x) const {
		solveLower(ret, x);
		solveLowerTranspose_inplace(ret);
	}

private : 
	void solveLower(VectorN<T>& ret, const VectorN<T> &x) const;
	void solveLowerTranspose_inplace(VectorN<T> &x) const;
	
public : 
	size_t						n;

	std::vector<T>				invdiag;	// reciprocals of diagonal elements

	std::vector<T>				value;		// values below the diagonal, listed column by column
	std::vector<size_t>			rowindex;	// a list of all row indices, for each column in turn
	std::vector<size_t>			colstart;	// where each column begins in rowindex (plus an extra entry at the end, of #nonzeros)

	std::vector<T>				adiag;		// just used in factorization: minimum "safe" diagonal entry allowed
};

template <class T> class PCGSolver {
public : 
	PCGSolver()									{ setParameters(1e-6, 300);							}
	PCGSolver(const SparseMatrix<T>& mat)		{ setParameters(1e-6, 300);	setSparseMatrix(mat);	}	

	void setSparseMatrix(const SparseMatrix<T>& mat) { fixed_matrix.set(mat); allocateMemory(); makePreconditioner_IC0(mat);	}	
	void allocateMemory() { 
		if(s.size() != fixed_matrix.n) { 
			size_t n = fixed_matrix.n; 
			s.resize(n); z.resize(n); r.resize(n); 
		} 
	}
	void makePreconditioner_IC0(const SparseMatrix<T>& mat) {
		ic_factor.makePreconditioner_IC0(mat);		
	}

	void setParameters(T _tol, size_t _maxIter, T _modified_incomplete_cholesky_parameter = 0.97f, T _min_diagonal_ratio = 0.25f) {
		tolerance = (_tol > 1.0e-30) ? _tol : 1.0e-30;
		maxIter = _maxIter;
		modified_incomplete_cholesky_parameter = _modified_incomplete_cholesky_parameter;
		min_diagonal_ratio = _min_diagonal_ratio;
	}

	bool cgSolve(VectorN<T> &result, const VectorN<T> &rhs, T &residual_out, size_t &iterations_out);
	bool pcgSolve(VectorN<T> &result, const VectorN<T> &rhs, T &residual_out, size_t &iterations_out);

protected:
	FixedSparseMatrix<T>		fixed_matrix;		// used within loop
	SparseColumnLowerFactor<T>	ic_factor;			// modified incomplete cholesky factor
	VectorN<T>					z, s, r;			// temporary vectors for PCG	

	// parameters
	T			tolerance;
	size_t		maxIter;
	T			modified_incomplete_cholesky_parameter;
	T			min_diagonal_ratio;
};

template<typename T>
SparseColumnLowerFactor<T>::SparseColumnLowerFactor(size_t _n=0) : n(_n), invdiag(_n), colstart(_n+1), adiag(_n) {}

template<typename T>
void SparseColumnLowerFactor<T>::clear() { resize(0); }

template<typename T>
void SparseColumnLowerFactor<T>::resize(size_t _n) {
	n = _n;
	invdiag.resize(n);
	value.clear();
	rowindex.clear();
	colstart.resize(n+1);
	adiag.resize(n);
}

template<typename T>
void SparseColumnLowerFactor<T>::makePreconditioner_IC0(const SparseMatrix<T>& matrix, T modification_parameter, T min_diagonal_ratio) {
	
	// first copy lower triangle of matrix into factor (Note: assuming A is symmetric of course!)
	resize(matrix.n);
	basis::vectorLibrary::zero(invdiag);		// important: eliminate old values from previous solves!
	basis::vectorLibrary::zero(adiag);

	for(size_t i=0; i<n; ++i) {
		colstart[i] = rowindex.size();
		for(size_t j=0; j<matrix.index[i].size(); ++j) {
			if(matrix.index[i][j] > i) {
				rowindex.push_back(matrix.index[i][j]);
				value.push_back(matrix.value[i][j]);
			}
			else if(matrix.index[i][j] == i) {
				invdiag[i] = adiag[i] = matrix.value[i][j];
			}
		}
	}
	colstart[matrix.n] = rowindex.size();
	
	// Now, Do the incomplete factorization (figure out numerical values)
	for(size_t k=0; k<matrix.n; ++k) {
		if(adiag[k] == 0)	// null row/column
			continue; 

		// figure out the final L(k,k) entry
		if(invdiag[k] < min_diagonal_ratio * adiag[k])
			invdiag[k] = 1 / sqrt(adiag[k]);			// drop to Gauss-Seidel here if the pivot looks dangerously small
		else
			invdiag[k] = 1 / sqrt(invdiag[k]);

		//invdiag[k] = 1.0 / sqrtf((float)::abs(invdiag[k])+0.001f);

		// finalize the k th column L(:,k)
		for(size_t p=colstart[k]; p<colstart[k+1]; ++p) {
			value[p] *= invdiag[k];
		}

		// incompletely eliminate L(:,k) from future columns, modifying diagonals
		for(size_t p=colstart[k]; p<colstart[k+1]; ++p) {
			size_t j = rowindex[p];		// work on column j
			invdiag[j] -= value[p] * value[p];
		}
	}
}

template<typename T>
void SparseColumnLowerFactor<T>::solveLower(VectorN<T>& ret, const VectorN<T> &x) const {
	assert(n == x.dim() && n == ret.dim());

	ret = x;
	for(size_t i=0; i<n; ++i) {
		ret[i] *= invdiag[i];
		for(size_t j = colstart[i]; j<colstart[i+1]; ++j) {
			ret[rowindex[j]] -= value[j] * ret[i];
		}
	}
}

template<typename T>
void SparseColumnLowerFactor<T>::solveLowerTranspose_inplace(VectorN<T> &x) const {
	assert(n == x.dim());

	size_t i = n;
	do {
		--i;
		for(size_t j = colstart[i]; j < colstart[i+1]; ++j) {
			x[i] -= value[j] * x[rowindex[j]];
		}
		x[i] *= invdiag[i];
	} while(i != 0);
}

template<typename T>
bool PCGSolver<T>::cgSolve(VectorN<T> &result, const VectorN<T> &rhs, T &residual_out, size_t &iterations_out) {

	size_t n = fixed_matrix.n;

	result.zero();
	r = rhs;

	residual_out = lInfNorm(r);
	if(residual_out == 0) {
		iterations_out = 0;
		return true;
	}

	T tol = tolerance * residual_out;

	s = r;
	T rho = dot(r, r);	
	
	size_t iteration;
	T alpha, beta, rho_new;

	for(iteration = 0; iteration < maxIter; ++iteration) {
		mul(z, fixed_matrix, s);			// z = A*s
		alpha = rho / dot(s, z);			// alpha = sigma / s.z
		result += s * alpha;				// p += alpha*s
		r -= z * alpha;						// r -= alpha*z	

		residual_out = lInfNorm(r);
		if(residual_out <= tolerance) {
			iterations_out = iteration+1;
			return true; 
		}

		rho_new = dot(r, r);				// sigma_new = r.r
		beta = rho_new / rho;				// beta = sigma_new / rho
		s *= beta; s += r;					// s = r + beta*s
		rho = rho_new;						// sigma = sigma_new
	}

	iterations_out = iteration;
	return false;
}

template<typename T>
bool PCGSolver<T>::pcgSolve(VectorN<T> &result, const VectorN<T> &rhs, T &residual_out, size_t &iterations_out) {

	size_t n = fixed_matrix.n;

	result.zero();
	r = rhs;

	residual_out = lInfNorm(r);
	if(residual_out == 0) {
		iterations_out = 0;
		return true;
	}

	//T tol = tolerance * residual_out;
	T tol = tolerance;

	ic_factor.applyPreconditioner(s, r);
	
	T rho = dot(s, r);
	if(rho == 0) {
		iterations_out = 0;
		return true;
	}
	
	size_t iteration;
	T alpha, beta, rho_new;

	for(iteration=0; iteration < maxIter; ++iteration) {
		mul(z, fixed_matrix, s);
		alpha = rho / dot(s, z);
		result += s * alpha;
		r -= z * alpha;

		residual_out = lInfNorm(r);
		if(residual_out <= tolerance) {
			iterations_out = iteration+1;
			return true; 
		}

		ic_factor.applyPreconditioner(z, r);
		rho_new = dot(z, r);
		beta = rho_new/rho;
		s *= beta; s += z;
		rho = rho_new;
	}

	iterations_out = iteration;
	return false;
}

};		// arithmetic

#endif
