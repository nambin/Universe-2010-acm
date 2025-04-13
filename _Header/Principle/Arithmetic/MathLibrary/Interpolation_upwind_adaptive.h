#pragma once

#include "../../../Basis/Random/RandomGen.h"
#include "Interpolation_upwind.h"

namespace arithmetic {

template<typename T>
const T cd2_adative(const T v[3], const T h[2]) {
	T dev[2], dev2;

	for(std::size_t i=0;i<2;++i) dev[i] = (v[i+1] - v[i]) / h[i];
	dev2 = (dev[1] - dev[0]) / (h[0] + h[1]);

	return dev[0] + h[0] * dev2;
}

template<typename T>
void eno2nd_adative(T& dx_m, T& dx_p, const T v[5], const T h[4]) {
	T dev[4], dev2[3];

	for(std::size_t i=0;i<4;++i) dev[i] = (v[i+1] - v[i]) / h[i];
	for(std::size_t i=0;i<3;++i) dev2[i] = (dev[i+1] - dev[i]) * T(2.0) / (h[i] + h[i+1]);

	dx_m = dev[1] + h[1] * T(0.5) * minmod(dev2[0], dev2[1]);
	dx_p = dev[2] - h[2] * T(0.5) * minmod(dev2[1], dev2[2]);
}

// Get high order 1st derivative through combining 1st upwind derivative
// for ith vertex (ith = 0, 1, ... , n)
template<typename T>
void combining_weight_1st_derivative(T* weight, const std::size_t ith, const std::size_t n, const T* h) {					
	// weight, h : size n array
	// ith = 0, 1, ... , n

	T* x = new T[n+1];
	T* w_x = new T[n+1];

	x[0] = 0;
	for(size_t i=0;i<n;++i) 
		x[i+1] = x[i] + h[i];		

	lagrange_polynomial_derivative(w_x, x[ith], n+1, x);

	weight[0] = - w_x[0] * h[0];
	for(size_t i=1;i<n;++i) 
		weight[i] = - (w_x[i] - weight[i-1] / h[i-1]) * h[i];

	delete[] x; delete[] w_x; 
}

template<typename T>
void weno_weight_from_eno(T weight[3], const std::size_t ith, const T h[5]) {		// ith = 2 or 3
	if(ith != 2 && ith != 3) {
		cout << "Error in weno_weight_from_eno(T weight[3], const std::size_t ith, const T h[5])" << endl;
		system("pause");
	}

	VectorN<T> weno(5);
	combining_weight_1st_derivative((T*)weno.getptr(), ith, 5, h);

	MatrixN<T> mat(5, 3); mat.zero();

	T * eno = new T[3];
	T * h_array = new T[3];
	for(std::size_t i=0;i<3;++i) {
		for(std::size_t k=0;k<3;++k)
			h_array[k] = h[i + k];

		combining_weight_1st_derivative(eno, ith - i, 3, h_array); 
		for(std::size_t k=0;k<3;++k)
			mat(i+k, i) = eno[k];
	}

	VectorN<T> weight_vector(3);
	SVD<T>(mat).solve(weight_vector, weno);

	for(size_t i=0;i<3;++i) 
		weight[i] = weight_vector[i];

	delete[] eno;
	delete[] h_array;
}

template<typename T>
void weno3rd_weight_from_eno(T weight[2], const std::size_t ith, const T h[3]) {		// ith = 1 or 2
	if(ith != 1 && ith != 2) {
		cout << "Error in weno3rd_weight_from_eno(T weight[2], const std::size_t ith, const T h[3])" << endl;
		system("pause");
	}

	VectorN<T> weno(3);
	combining_weight_1st_derivative((T*)weno.getptr(), ith, 3, h);

	MatrixN<T> mat(3, 2); mat.zero();

	T * eno = new T[2];
	T * h_array = new T[2];
	for(std::size_t i=0;i<2;++i) {
		for(std::size_t k=0;k<2;++k)
			h_array[k] = h[i + k];

		combining_weight_1st_derivative(eno, ith - i, 2, h_array); 
		for(std::size_t k=0;k<2;++k)
			mat(i+k, i) = eno[k];
	}

	VectorN<T> weight_vector(2);
	SVD<T>(mat).solve(weight_vector, weno);

	for(size_t i=0;i<2;++i) 
		weight[i] = weight_vector[i];

	delete[] eno;
	delete[] h_array;
}

};			// arithmetic

#include "Interpolation_weno3rd.h"
#include "Interpolation_weno5th.h"