#ifndef _ARITHMETIC_INTERPOLATION_H_
#define _ARITHMETIC_INTERPOLATION_H_

#include <iostream>
#include <cmath>

namespace arithmetic {

template <typename T, typename S> 
inline const T lerp(const T &c0, const T &c1, S weight0)	{
	return c0 * weight0 + c1 * (S(1.0) - weight0);;
}

template <typename T, typename S> 
inline void lerp(T &interp, const T &c0, const T &c1, S weight0)	{
	interp = c0 * weight0 + c1 * (S(1.0) - weight0);;
}

template <typename T, typename S> 
inline const T lerp(const T &c0, const T &c1, S weight0, S weight1)	{
	return (c0 * weight0 + c1 * weight1) / (weight0 + weight1);
}

template <typename T, typename S> 
inline void lerp(T &interp, const T &c0, const T &c1, S weight0, S weight1)	{
	interp = (c0 * weight0 + c1 * weight1) / (weight0 + weight1);
}

template <typename T>
void lagrange_polynomial(T* weight, const T x, std::size_t n, const T* pos_sample) {		
	// weight, pos_sample : size n array

	for(std::size_t i=0;i<n;++i) {
		weight[i] = 1;
		for(std::size_t k=0;k<n;++k) {
			if(i == k)
				continue;

			T nomi = (		x		- pos_sample[k]);
			T deno = (pos_sample[i] - pos_sample[k]);

			weight[i] *= ( nomi / deno );
		}
	}
}

template <typename T>
void lagrange_polynomial_derivative(T* weight, const T x, std::size_t n, const T* pos_sample) {		
	// weight, pos_sample : size n array
	
	for(std::size_t i=0;i<n;++i) {
		weight[i] = 0;

		T denom = 1;
		for(std::size_t k=0;k<n;++k) {
			if(i == k)
				continue;

			denom *= (pos_sample[i] - pos_sample[k]);			
		}

		T sum = 0;
		for(std::size_t k=0;k<n;++k) {
			if(i == k)
				continue;

			T _mul = 1;
			for(std::size_t k_=0;k_<n;++k_) {
				if(k_ == k || k_ == i)
					continue;

				_mul *= (x - pos_sample[k_]);
			}

			sum += _mul;			
		}

		weight[i] += ( sum / denom );
	}
}

//template <typename T, typename S> 
//inline void cubicInterpolation( T f0, T f1, T d0, T d1, S x, S h, T &f, T &d )	{
//
//	T delta = f1 - f0;
//
//	d0 *= h;
//	d1 *= h;
//
//	// Make monotonic
//	//if ( delta == 0.0 ) {
//	//	d0 = 0.0; d1 = 0.0;
//	//} else if ( delta > 0.0 ) {
//	//	if ( d0 < 0.0 ) d0 = 0.0;
//	//	else if ( d0 > 3.0*delta ) 
//	//		d0 = 3.0*delta;
//	//	if ( d1 < 0.0 )	d1 = 0.0;
//	//	else if ( d1 > 3.0*delta )
//	//		d1 = 3.0*delta;
//	//} else { // delta < 0.0
//	//	if ( d0 > 0.0 ) d0 = 0.0;
//	//	else if ( d0 < 3.0*delta ) 
//	//		d0 = 3.0*delta;
//	//	if ( d1 > 0.0 ) d1 = 0.0;
//	//	else if ( d1 < 3.0*delta ) 
//	//		d1 = 3.0*delta;
//	//}
//
//	T a = (d0 + d1) - 2.0f*delta;
//	T b = 3.0f*delta - (2.0f*d0 + d1);
//
//	f = ((a*x + b)*x + d0)*x + f0;
//	d = ((3.0f*a*x + 2.0f*b)*x + d0)/h;
//}

};		

#include "Interpolation_upwind.h"
#include "Interpolation_upwind_adaptive.h"
#include "Interpolation_doyub.h"

#endif