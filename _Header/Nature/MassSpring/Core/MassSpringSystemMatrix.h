#ifndef _MASS_SPRING_SYSTEM_MATRIX_H_
#define _MASS_SPRING_SYSTEM_MATRIX_H_

#include <vector>
#include "../../../Basis/Basis/Array.h"
#include "../../../Principle/Arithmetic/Matrix/Matrix3.h"
#include "../../../Principle/Arithmetic/Vector/Vector3.h"
#include "../Core/SpringLibrary.h"

namespace springspace {

using std::vector;
using basis::Array;
using arithmetic::Vector3;
using arithmetic::Matrix3;

template<typename T>
class MassSpringSystemMatrix {
public : 
	MassSpringSystemMatrix();

	void clear();
	void zero();
	size_t numParticles() const { return n; }

	template <typename T> friend void mul(Vec3Array<T> &ret, const MassSpringSystemMatrix<T> &mat, const Vec3Array<T> &v);

	void make_MatrixIndexMapping(Array<std::pair<size_t, size_t>>& mapping, const Array<Spring<T>>& springs) const;

public : 
	size_t					n;			// number of particles

	vector<Matrix3<T>>		diag;

	vector<Matrix3<T>>		value;
	vector<size_t>			colindex;	
	vector<size_t>			rowstart;	
};

template<typename T>
MassSpringSystemMatrix<T>::MassSpringSystemMatrix() { clear(); }

template<typename T>
void MassSpringSystemMatrix<T>::clear() {
	n = 0;
	diag.clear();
	value.clear();
	colindex.clear();
	rowstart.clear();
}

template<typename T>
void MassSpringSystemMatrix<T>::zero() {
	for(size_t i=0;i<diag.size();++i)	diag[i].zero();
	for(size_t i=0;i<value.size();++i)	value[i].zero();
}

template <typename T> 
void mul(Vec3Array<T> &ret, const MassSpringSystemMatrix<T> &mat, const Vec3Array<T> &v) {
	assert(ret.size() == v.size() && ret.size() == mat.numParticles());

	for(size_t i=0;i<mat.n;++i) {
		ret[i] = mat.diag[i] * v[i];
		for(size_t j=mat.rowstart[i];j<mat.rowstart[i+1];++j)
			ret[i] += mat.value[j] * v[mat.colindex[j]];
	}
}

template <typename T> 
void MassSpringSystemMatrix<T>::make_MatrixIndexMapping(Array<std::pair<size_t, size_t>>& mapping, const Array<Spring<T>>& springs) const {

	size_t idx0, idx1;

	mapping.resize(springs.size());
	for(size_t i=0;i<springs.size();++i) {
		idx0 = springs[i].idx0;
		idx1 = springs[i].idx1;

		for(size_t k=rowstart[idx0];k<rowstart[idx0+1];++k) {
			if(idx1 == colindex[k]) {
				mapping[i].first = k;
				break;
			}
		}
		for(size_t k=rowstart[idx1];k<rowstart[idx1+1];++k) {
			if(idx0 == colindex[k]) {
				mapping[i].second = k;
				break;
			}
		}
	}
}

};		// springspace

#endif