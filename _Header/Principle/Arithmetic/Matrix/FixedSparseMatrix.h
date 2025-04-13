#ifndef _ARITHMETIC_FIXEDSPARSEMATRIX_H_
#define _ARITHMETIC_FIXEDSPARSEMATRIX_H_

#include "SparseMatrix.h"

namespace arithmetic {

template<class T> class FixedSparseMatrix {
public : 
	typedef std::size_t			index_type;
	typedef index_type*			index_pointer ;
	typedef const index_type*	index_const_pointer ;
	typedef index_type&			index_reference ;
	typedef const index_type&	index_const_reference ;

	explicit FixedSparseMatrix(index_type n_=0);
	explicit FixedSparseMatrix(const SparseMatrix<T>&);
	void set(const SparseMatrix<T>&);

	void clear();
	void resize(index_type);
	
	index_type dim() const { return n; }

	template <typename T> friend void mul(VectorN<T> &ret, const FixedSparseMatrix<T> &m, const VectorN<T> &v);

public : 
	index_type						n;			// dimension

	std::vector<T>					value;		// nonzero values row by row
	std::vector<index_type>			colindex;	// corresponding column indices
	std::vector<index_type>			rowstart;	// where each row starts in value and colindex (and last entry is one past the end, the number of nonzeros)
};

template<typename T>
FixedSparseMatrix<T>::FixedSparseMatrix(index_type _n) { resize(_n); }

template<typename T>
FixedSparseMatrix<T>::FixedSparseMatrix(const SparseMatrix<T>& mat) { set(mat); }

template<typename T>
void FixedSparseMatrix<T>::set(const SparseMatrix<T>& mat) {
	resize(mat.n);
	rowstart[0] = 0;
	for(index_type i=0; i<n; ++i) {
		rowstart[i+1] = rowstart[i] + mat.index[i].size();
	}

	value.resize(rowstart[n]);
	colindex.resize(rowstart[n]);
	index_type j=0;
	for(index_type i=0; i<n; ++i) {
		for(index_type k=0; k<mat.index[i].size(); ++k) {
			value[j] = mat.value[i][k];
			colindex[j] = mat.index[i][k];
			++j;
		}
	}
}

template<typename T>
void FixedSparseMatrix<T>::clear() {
	resize(0);
}

template<typename T>
void FixedSparseMatrix<T>::resize(index_type n_) {
	n = n_;
	value.clear();
	colindex.clear();
	rowstart.resize(n+1);
}

template <typename T> 
void mul(VectorN<T> &ret, const FixedSparseMatrix<T> &mat, const VectorN<T> &v) {
	assert(mat.dim() == v.dim() && ret.dim() == v.dim());

	for(FixedSparseMatrix<T>::index_type i=0;i<mat.n;++i) {
		ret[i] = 0;
		for(FixedSparseMatrix<T>::index_type j=mat.rowstart[i];j<mat.rowstart[i+1];++j)
			ret[i] += mat.value[j] * v[mat.colindex[j]];

		//for(index_type j=0;j<5;++j) ret[i] += mat.value[j];
		//ret[i] += (mat.value[0]+mat.value[1]+mat.value[2]+mat.value[3]+mat.value[4]);
	}
}

};

#endif
