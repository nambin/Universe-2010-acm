#ifndef _ARITHMETIC_PCA_H_
#define _ARITHMETIC_PCA_H_

#include <vector>
#include "../../Principle.h"
#include "SymmEig.h"

namespace arithmetic {

#ifdef _DEBUG
#define __ASSERT(c) assert(c)
#else
#define __ASSERT(c) {}
#endif

template<typename T> class PCA {
public : 
	typedef typename MatrixN<T>::index_type index_type;

	//////////////////////////////////////////////////////////////////////////
	// Constructors
	PCA(index_type dim);
	void init(index_type dim);

	//////////////////////////////////////////////////////////////////////////
	// Manipulation
	void push(const VectorN<T> &v);
	void doPCA();

	//////////////////////////////////////////////////////////////////////////
	// Getter
	inline index_type getDim() const		{ return _dim;		}
	inline index_type getNumData() const	{ return _numData;	}
	inline index_type getNumBasis() const	{ return _numBasis; }	

	void getAverage(VectorN<T> &avg) const;
	inline void getBasis(index_type ith, VectorN<T> &basis) const;
	inline const T getVariance(index_type ith) const;
	inline const T getDeviation(index_type ith) const;
	
	inline void getCoefficient(index_type iData, VectorN<T> &coeff) const;
	void getCoefficient(const VectorN<T> &v, VectorN<T> &coeff) const;

	//////////////////////////////////////////////////////////////////////////
	// IO
	template <typename T> friend std::ostream& operator<<(std::ostream &os, const PCA<T> &pca);	

private :
	//////////////////////////////////////////////////////////////////////////
	// IO
	void calcAverage();

private : 
	index_type _dim, _numData, _numBasis;
	std::vector<VectorN<T>> data;

	VectorN<T> _average;
	std::vector<VectorN<T>> _basis;
	std::vector<T>			_variance;
};

//////////////////////////////////////////////////////////////////////////
// Constructors
template<typename T> PCA<T>::PCA(index_type dim) {
	init(dim);
}

template<typename T> void PCA<T>::init(index_type dim) {
	__ASSERT(dim > 0);
	_dim = dim;
	_numData = 0;
	_numBasis = 0;
	data.clear();
	_average.resize(_dim);
	_basis.clear();
	_variance.clear();
}

//////////////////////////////////////////////////////////////////////////
// Manipulation
template<typename T> void PCA<T>::push(const VectorN<T> &v) {
	__ASSERT(_dim == v.dim());

	data.push_back(v);
	_numData = static_cast<index_type>(data.size());
}

template<typename T> void PCA<T>::doPCA() {

	VectorN<T> avg(_dim), tmp(_dim);
	calcAverage();
	getAverage(avg);

	_numData = static_cast<index_type>(data.size());

	MatrixN<T> matData(_dim, _numData);
	for(index_type i=0;i<_numData;++i) {
		tmp = data[i] - avg;
		matData.setCol(i, tmp);
	}															// matData = data - average

	if(_dim <= _numData) {
		MatrixN<T> mat(_dim, _dim);
		for(index_type i=0;i<_dim;++i) 
			for(index_type j=0;j<_dim;++j) 
				for(index_type k=0;k<_numData;++k) 
					mat(i,j) += matData(i,k) * matData(j,k);		// mat = matData * matData^T

		SymmEig<T> symmeig(mat);									// Perform SVD

		/////////////////////////
		// Set basis, variance //
		/////////////////////////
		_basis.resize(_dim);		for(index_type i=0;i<_dim;++i) _basis[i].resize(_dim);
		_variance.resize(_dim);

		for(index_type i=0;i<_dim;++i) {
			symmeig.z.getCol(i, _basis[i]);
			_variance[i] = symmeig.d[i] / static_cast<T>(_numData-1);
		}
	}
	else {
		if(_numData == 1) 
			return;

		MatrixN<T> mat(_numData, _numData);
		for(index_type i=0;i<_numData;++i) 
			for(index_type j=0;j<_numData;++j) 
				for(index_type k=0;k<_dim;++k) 
					mat(i,j) += matData(k,i) * matData(k,j);		// mat = matData^T * matData

		SymmEig<T> symmeig(mat);									// Perform SVD

		/////////////////////////
		// Set basis, variance //
		/////////////////////////
		_basis.resize(_numData);		for(index_type i=0;i<_numData;++i) _basis[i].resize(_dim);
		_variance.resize(_numData);

		VectorN<T> _eigenVector(_numData);
		for(index_type i=0;i<_numData;++i) {
			symmeig.z.getCol(i, _eigenVector);
			for(index_type j=0;j<_numData;++j) 
				_basis[i] += _eigenVector[j] * data[j];
			normalize(_basis[i]);
			_variance[i] = symmeig.d[i] / static_cast<T>(_numData-1);
		}
	}	

	_numBasis = MIN(_dim, _numData - 1);
}

template<typename T> void PCA<T>::calcAverage() {

	_numData = static_cast<index_type>(data.size());
	for(index_type i=0;i<_numData;++i) 
		_average += data[i];

	_average /= static_cast<T>(_numData);
}

//////////////////////////////////////////////////////////////////////////
// Getter
template<typename T> void PCA<T>::getAverage(VectorN<T> &avg) const {
	avg = _average;
}

template<typename T> inline void PCA<T>::getBasis(index_type ith, VectorN<T> &basis) const {
	basis = _basis[ith];
}

template<typename T> inline const T PCA<T>::getVariance(index_type ith) const {
	return _variance[ith];
}

template<typename T> inline const T PCA<T>::getDeviation(index_type ith) const {
	if(_variance[ith] < 0)
		return 0;
	else 
		return sqrt(_variance[ith]);
}

template<typename T> inline void PCA<T>::getCoefficient(index_type iData, VectorN<T> &coeff) const {
	getCoefficient(data[iData], coeff);
}

template<typename T> void PCA<T>::getCoefficient(const VectorN<T> &v, VectorN<T> &coeff) const {
	index_type dim = coeff.dim();
	VectorN<T> _diff(dim);
	_diff = v - _average;

	for(index_type i=0;i<dim;++i)
		coeff[i] = dot(_diff, _basis[i]);
}

//////////////////////////////////////////////////////////////////////////
// IO
template <typename T> std::ostream& operator<<(std::ostream &os, const PCA<T> &pca) {

	typename PCA<T>::index_type numComp = 5;

	os << "# Data Dimension : " << pca._dim << std::endl;
	os << "# Data Average" << std::endl;

	if(pca._dim < numComp) os << pca._average << std::endl;
	else {
		os << "[";
		for(typename PCA<T>::index_type i=0;i<numComp;++i) 
			os << pca._average[i] << " ";
		os << "...]" << std::endl;
	}

	os << "# Number of Basis : " << pca._numBasis << std::endl;
	os << "# Basis _ Variance _ Standard Deviation" << std::endl;
	for(typename PCA<T>::index_type i=0;i<pca._numBasis;++i) {
		os << i << "th _ ";

		if(pca._dim < numComp) os << pca._basis[i] << std::endl;
		else {
			os << "[";
			for(typename PCA<T>::index_type k=0;k<numComp;++k) 
				os << pca._basis[i][k] << " ";
			os << "...] _ ";
		}

		os << pca.getVariance(i) << " _ " << pca.getDeviation(i) << std::endl;
	}

	return os;
}

void PRINCIPLEDLL testPCA();

#undef __ASSERT

};		// namespace arithmetic

#endif

