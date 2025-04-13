#ifndef _MASS_SPRING_VEC3_ARRAY_H_
#define _MASS_SPRING_VEC3_ARRAY_H_

#include "../../../Basis/Basis/Array.h"
#include "../../../Principle/Arithmetic/Vector/Vector3.h"

namespace springspace {

using basis::Array;
using arithmetic::Vector3;

template<typename T>
class Vec3Array : public Array<Vector3<T>> {
public : 
	Vec3Array();
	Vec3Array(const Vec3Array &src);
	explicit Vec3Array(index_type size);
	~Vec3Array();
};

//////////////////////////////////////////////////////////////////////////
// Constructors
template <typename T> Vec3Array<T>::Vec3Array() : Array<Vector3<T>>() {}
template <typename T> Vec3Array<T>::Vec3Array(const Vec3Array &src) : Array<Vector3<T>>(src) {}
template <typename T> Vec3Array<T>::Vec3Array(index_type size) : Array<Vector3<T>>(size) {}
template <typename T> Vec3Array<T>::~Vec3Array() {}

template<typename T> const T lInfNorm(const Vec3Array<T>& x) {
	T max = -1;	
	for(Vec3Array<T>::index_type i=0;i<x.size();++i) {
		max = (max > ::abs(x[i][0])) ? max : ::abs(x[i][0]);
		max = (max > ::abs(x[i][1])) ? max : ::abs(x[i][1]);
		max = (max > ::abs(x[i][2])) ? max : ::abs(x[i][2]);
	}
	return max;
}

template<typename T> const T dot(const Vec3Array<T>& x, const Vec3Array<T>& y) {
	T _dot = 0;	
	for(Vec3Array<T>::index_type i=0;i<x.size();++i) {
		_dot += (x[i][0] * y[i][0]);
		_dot += (x[i][1] * y[i][1]); 
		_dot += (x[i][2] * y[i][2]);
	}
	return _dot;
}

template<typename T> void increment(Vec3Array<T>& x, T scale, const Vec3Array<T>& y) {
	for(Vec3Array<T>::index_type i=0;i<x.size();++i) 
		x[i] += (scale * y[i]);
}

template<typename T> void printMaxVector3(Vec3Array<T>& vec) {
	T sqrmag = -1;
	size_t idx;
	for(size_t i=0;i<vec.size();++i) {
		if(sqrmag < arithmetic::sqrmag(vec[i])) {
			sqrmag = arithmetic::sqrmag(vec[i]);
			idx = i;
		}
	}
	std::cout << idx << " : " << vec[idx] << std::endl;
}


};			// springspace

#endif