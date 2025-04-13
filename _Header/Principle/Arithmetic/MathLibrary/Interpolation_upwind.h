#pragma once

#include "../Matrix/MatrixN.h"
#include "SVD.h"
#include "Functions.h"
#include "Interpolation_upwind_adaptive.h"

namespace arithmetic {

template<typename T>
const T minmod(const T x, const T y) {
	if( x * y < 0 )								return 0;
	else if( ::abs(x) < ::abs(y) )				return x;
	else										return y;
}

template<typename T>
const T eno2nd(const T v[5], const T h, bool is_velocity_positive) {	// i-2, i-1, i, i+1, i+2

	T hInv = T(1.0/h), hInvhalf = hInv*T(0.5);
	T q1, c;

	int i = 2;
	int k = i;								// if vel < 0, then k = i
	if(is_velocity_positive) k = i-1;		// if vel > 0, then k = i-1

	T v_1_dev[3];															// k-1_2, k+1_2, k+3_2
	for(int m=0;m<3;++m) v_1_dev[m] = (v[k+m] - v[k-1+m])*hInv;

	q1 = v_1_dev[1];						// First Order _ D1(k+1_2)

	T v_2_dev[2];															// k, k+1
	for(int m=0;m<2;++m) v_2_dev[m] = (v_1_dev[m+1] - v_1_dev[m])*hInvhalf;

	// minmod limiter
	c = minmod(v_2_dev[0], v_2_dev[1]);	

	return q1 + c * T(2*(i-k) - 1) * h;
}

template<typename T>
const T eno(const T v[7], const T h, bool is_velocity_positive) {	// i-3, i-2, i-1, i, i+1, i+2, i+3

	T hInv = T(1.0/h), hInvhalf = hInv*T(0.5), hInvthird = T(hInv/3.0);
	T q1, c, cstar;

	int i = 3;
	int kstar, k = i;						// if vel < 0, then k = i
	if(is_velocity_positive) k = i-1;		// if vel > 0, then k = i-1

	T v_1_dev[5];															// k-3_2 , k-1_2, k+1_2, k+3_2, k+5_2
	for(int m=0;m<5;++m) v_1_dev[m] = (v[k-1+m] - v[k-2+m])*hInv;

	q1 = v_1_dev[2];				// First Order _ D1(k+1_2)

	T v_2_dev[4];															// k-1, k, k+1, k+2
	for(int m=0;m<4;++m) v_2_dev[m] = (v_1_dev[m+1] - v_1_dev[m])*hInvhalf;

	if(::abs(v_2_dev[1]) < ::abs(v_2_dev[2]))	kstar = k-1;		// if |D2k| < |D2k+1|, then kstar = k-1
	else										kstar = k;			// if |D2k| > |D2k+1|, then kstar = k

	int kstar_local = kstar - (k - 1);
	c = v_2_dev[kstar_local+1];		// Second Order	_ D2(kstar+1)

	T v_3_dev[2];												// kstar+1_2, kstar+3_2
	v_3_dev[0] = (v_2_dev[kstar_local+1] - v_2_dev[kstar_local]) * hInvthird;
	v_3_dev[1] = (v_2_dev[kstar_local+2] - v_2_dev[kstar_local+1]) * hInvthird;

	// Third Order
	if(::abs(v_3_dev[0]) < ::abs(v_3_dev[1]))	cstar = v_3_dev[0];
	else										cstar = v_3_dev[1];

	return q1 + c * T(2*(i-k) - 1) * h + cstar * (3*SQ(i-kstar) - 6*(i-kstar) + 2) * SQ(h);
}

template<typename T>
const T weno(const T v[7], const T h, bool is_velocity_positive, const T eps = 1.0e-8) {	// i-3, i-2, i-1, i, i+1, i+2, i+3

	static const T c_1_3 = T(1.0/3.0), c_1_4 = T(0.25), c_1_6 = T(1.0/6.0), c_5_6 = T(5.0/6.0), c_7_6 = T(7.0/6.0), c_11_6 = T(11.0/6.0);
	static const T c_13_12 = T(13.0/12.0);

	T hInv = T(1.0/h);	
	T vdev[5];

	if(is_velocity_positive) { for(int m=0;m<5;++m) vdev[m] = (v[m+1] - v[m  ]) * hInv;	}
	else					 { for(int m=0;m<5;++m) vdev[m] = (v[6-m] - v[5-m]) * hInv; }

	//T max_vsq = 0; for(int m=0;m<5;++m) if(SQ(vdev[m]) > max_vsq) max_vsq = SQ(vdev[m]);
	//T eps = 1.0e-6 * max_vsq + 1.0e-8;

	T phix1 =   vdev[0] * c_1_3  - vdev[1] * c_7_6 + vdev[2] * c_11_6;
	T phix2 = - vdev[1] * c_1_6  + vdev[2] * c_5_6 + vdev[3] * c_1_3;
	T phix3 =   vdev[2] * c_1_3  + vdev[3] * c_5_6 - vdev[4] * c_1_6;

	T s1 = c_13_12 * SQ(vdev[0] - 2*vdev[1] + vdev[2]) + c_1_4 * SQ(vdev[0] - 4*vdev[1] + 3*vdev[2]);
	T s2 = c_13_12 * SQ(vdev[1] - 2*vdev[2] + vdev[3]) + c_1_4 * SQ(vdev[1] - vdev[3]);
	T s3 = c_13_12 * SQ(vdev[2] - 2*vdev[3] + vdev[4]) + c_1_4 * SQ(3*vdev[2] - 4*vdev[3] + vdev[4]);

	T alpha1 = T(0.1 / SQ(s1 + eps));
	T alpha2 = T(0.6 / SQ(s2 + eps));
	T alpha3 = T(0.3 / SQ(s3 + eps));

	T sum = alpha1 + alpha2 + alpha3;

	return (alpha1 * phix1 + alpha2 * phix2 + alpha3 * phix3) / sum;
}

template<typename T>
const T weno3rd(const T v[5], const T h, bool is_velocity_positive, const T eps = 1.0e-8) {	// i-2, i-1, i, i+1, i+2

	static const T c_1_3 = T(1.0/3.0), c_2_3 = T(2.0/3.0);

	T hInv = T(1.0/h);	
	T vdev[3];

	if(is_velocity_positive) { for(int m=0;m<3;++m) vdev[m] = (v[m+1] - v[m  ]) * hInv;	}
	else					 { for(int m=0;m<3;++m) vdev[m] = (v[4-m] - v[3-m]) * hInv; }

	T phix1 = vdev[0] * T(-0.5) + vdev[1] * T(1.5);
	T phix2 = vdev[1] * T(+0.5) + vdev[2] * T(0.5);

	T s1 = SQ(vdev[0] - vdev[1]);
	T s2 = SQ(vdev[1] - vdev[2]);

	T alpha1 = T(c_1_3 / SQ(s1 + eps));
	T alpha2 = T(c_2_3 / SQ(s2 + eps));	

	T sum = alpha1 + alpha2;

	return (alpha1 * phix1 + alpha2 * phix2) / sum;
}

template<typename T>
void combining_weight_1st_derivative_uniform(T* weight, const std::size_t ith, const std::size_t n, const T h) {			
	// weight : size n array
	// ith = 0, 1, ... , n

	T* h_array = new T[n];
	for(size_t i=0;i<n;++i) 
		h_array[i] = h;

	combining_weight_1st_derivative(weight, ith, n, h_array);

	delete[] h_array;
}

template<typename T>
void weno_weight_from_eno_uniform(T weight[3], const std::size_t ith, const T h) {		// ith = 2 or 3
	T h_array[5]; for(int i=0;i<5;++i) h_array[i] = h;
	weno_weight_from_eno(weight, ith, h_array);
}

template<typename T>
void weno3rd_weight_from_eno_uniform(T weight[2], const std::size_t ith, const T h) {		// ith = 1 or 2
	T h_array[3]; for(int i=0;i<3;++i) h_array[i] = h;
	weno3rd_weight_from_eno(weight, ith, h_array);
}

};		// namespace arithmetic