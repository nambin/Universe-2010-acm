#ifndef _ARITHMETIC_SPARSEMATRIX_H_
#define _ARITHMETIC_SPARSEMATRIX_H_

#include <iostream>
#include "../../../Basis/vectorLibrary/vectorLibrary.h"
#include "../Vector/VectorN.h"

namespace arithmetic {

template<class T> class SparseMatrix { 
public : 
	typedef std::size_t			index_type;
	typedef index_type*			index_pointer ;
	typedef const index_type*	index_const_pointer ;
	typedef index_type&			index_reference ;
	typedef const index_type&	index_const_reference ;

	explicit SparseMatrix(index_type n_=0, index_type expected_nonzeros_per_row=7);	
	void clear();
	void zero();
	void resize(index_type);

	index_type dim() const { return n; }
	const T operator()(index_type i, index_type j) const;

	void set_element(index_type i, index_type j, T val);
	void add_element(index_type i, index_type j, T val);

	// Assumes indices are already sorted	
	void add_sparse_row(index_type i, const std::vector<index_type> &idxs, const std::vector<T> &vals); 

	// Assumes matrix has symmetric structure - so the indices in row i tell us which columns to delete i from
	void symmetric_remove_row_and_column(index_type i);

	template <typename T> friend void mul(VectorN<T> &ret, const SparseMatrix<T> &m, const VectorN<T> &v);

public : 
	index_type								n;		// dimension

	std::vector<std::vector<index_type>>	index;	// for each row, a list of all column indices (sorted)
	std::vector<std::vector<T>>				value;	// values corresponding to index
};

template<class T> 
SparseMatrix<T>::SparseMatrix(index_type n_, index_type expected_nonzeros_per_row) : n(n_), index(n_), value(n_) {
	for(index_type i=0;i<n;++i) {
		index[i].reserve(expected_nonzeros_per_row);
		value[i].reserve(expected_nonzeros_per_row);
	}
}

template<class T> 
void SparseMatrix<T>::clear() {
	n = 0;
	index.clear();
	value.clear();
}

template<class T> 
void SparseMatrix<T>::zero(void) {
	for(index_type i=0; i<n; ++i) {
		index[i].resize(0);
		value[i].resize(0);
	}
}

template<class T> 
void SparseMatrix<T>::resize(index_type n_) {
	n = n_;
	index.resize(n);
	value.resize(n);
}

template<class T> 
const T SparseMatrix<T>::operator()(index_type i, index_type j) const {
	for(index_type k=0; k<index[i].size(); ++k){
		if(index[i][k]==j)		return value[i][k];
		else if(index[i][k]>j)	return 0;
	}
	return 0;
}

template<class T> 
void SparseMatrix<T>::set_element(index_type i, index_type j, T val) {
	for(index_type k=0;k<index[i].size();++k) {
		if(index[i][k] == j) {
			value[i][k] = val;
			return;
		}
		else if(index[i][k] > j) {
			basis::vectorLibrary::insert(index[i], k, j);
			basis::vectorLibrary::insert(value[i], k, val);
			return;
		}
	}
	index[i].push_back(j);
	value[i].push_back(val);
}

template<class T> 
void SparseMatrix<T>::add_element(index_type i, index_type j, T val) {	
	for(index_type k=0; k<index[i].size(); ++k) {
		if(index[i][k] == j) {
			value[i][k] += val;
			return;
		}
		else if(index[i][k] > j) {
			basis::vectorLibrary::insert(index[i], k, j);
			basis::vectorLibrary::insert(value[i], k, val);
			return;
		}
	}
	index[i].push_back(j);
	value[i].push_back(val);
}

// Assumes indices are already sorted
template<class T> 
void SparseMatrix<T>::add_sparse_row(index_type i, const std::vector<index_type> &indices, const std::vector<T> &values) {
	index_type j=0, k=0;
	while(j<indices.size() && k<index[i].size()) {
		if(index[i][k] < indices[j]) {
			++k;
		}
		else if(index[i][k] > indices[j]) {
			basis::vectorLibrary::insert(index[i], k, indices[j]);
			basis::vectorLibrary::insert(value[i], k, values[j]);
			++j;
		}
		else {
			value[i][k] += values[j];
			++j;
			//++k;
		}
	}
	for(;j<indices.size(); ++j) {
		index[i].push_back(indices[j]);
		value[i].push_back(values[j]);
	}
}

// assumes matrix has symmetric structure - so the indices in row i tell us which columns to delete i from
template<class T> 
void SparseMatrix<T>::symmetric_remove_row_and_column(index_type i) {
	for(index_type a=0; a<index[i].size(); ++a) {
		index_type j = index[i][a];  
		for(index_type b=0; b<index[j].size(); ++b) {
			if(index[j][b] == i) {
				basis::vectorLibrary::erase(index[j], b);
				basis::vectorLibrary::erase(value[j], b);
				break;
			}
		}
	}
	index[i].resize(0);
	value[i].resize(0);
}

template <typename T> 
void mul(VectorN<T> &ret, const SparseMatrix<T> &m, const VectorN<T> &v) {
	assert(m.dim() == v.dim() && ret.dim() == v.dim());

	if(&ret == &v) {
		VectorN<T> _ret(ret.dim());

		mul(_ret, m, v);
		ret = _ret;

		return;
	}
	else {
		for(SparseMatrix<T>::index_type i=0; i<m.n; ++i) {
			ret[i] = 0;
			for(SparseMatrix<T>::index_type j=0; j<m.index[i].size(); ++j) {
				ret[i] += m.value[i][j] * v[m.index[i][j]];
			}
		}
	}	
}

};		// arithmetic

#endif
