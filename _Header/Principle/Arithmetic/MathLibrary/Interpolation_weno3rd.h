#pragma once

#include "Interpolation_upwind.h"

namespace arithmetic {

template<typename T>
class WENO3rd_Adaptive {
public : 
	WENO3rd_Adaptive() {}
	WENO3rd_Adaptive(const T h) { set(h); }
	WENO3rd_Adaptive(const T h[4]) { set(h); }

	void set(const T h) { const T _h[4] = {h,h,h,h}; pre_compute_weights(_h); }
	void set(const T h[4]) { pre_compute_weights(h); }

	void dx(T& dx_m, T& dx_p, const T v[5]) const;	
	void print_weights(const bool c = true, const bool s = true, const bool w = true) const;

	static void validate_weights_uniform_weno(const double h = 0.5);
	static void validate_weights_except_smoothness(const double h[4]);

	template<typename FUNC> static void validate_values_uniform_weno(const double h, const FUNC&, const double a, const double b, std::size_t n_sample = 10);

private :
	void pre_compute_weights(const T h[4]);

	static void calc_combining_weight_1st_derivative(T c[2][2], const T h[3], const bool is_positive_vel);
	static void calc_smoothness_weight(T b[2][2], const T h[4]);
	static void calc_weno_weight(T d[2], const T h[3], const bool is_positive_vel);

	const T dx(const T dev[3], const T c[2][2], const T s[2][2], const T d[2]) const;
	void print_weights(const T c[2][2], const T s[2][2], const T d[2], const bool c_w, const bool s_w, const bool w_w) const;

private : 
	T h[4], hInv[4];

	T c_m[2][2], c_p[2][2];

	T smooth_m[2][2];
	T smooth_p[2][2];	

	T d_m[2], d_p[2];
};

template<typename T>
void WENO3rd_Adaptive<T>::dx(T& dx_m, T& dx_p, const T v[5]) const {
	T dev[4];
	for(std::size_t i=0;i<4;++i) dev[i] = (v[i+1] - v[i]) * hInv[i];

	dx_m = dx(&dev[0], c_m, smooth_m, d_m);
	dx_p = dx(&dev[1], c_p, smooth_p, d_p);	
}

template<typename T> 
const T WENO3rd_Adaptive<T>::dx(const T dev[3], const T c[2][2], const T s[2][2], const T d[2]) const {
	
	static const T eps = 1.0e-8;

	T _dx[2] = {0,0}, _smooth[2] = {0,0}, _weight[2] = {0,0};	
	T s0;

	for(std::size_t i=0;i<2;++i) 
		for(std::size_t j=0;j<2;++j) 
			_dx[i] += (dev[1-i+j] * c[i][j]);

	for(std::size_t i=0;i<2;++i) {
		s0 = dev[1-i] * s[i][0] + dev[2-i] * s[i][1];			
		_smooth[i] = (s0 * s0);
	}			

	_weight[0] = d[0] / SQ(eps + _smooth[0]);
	_weight[1] = d[1] / SQ(eps + _smooth[1]);

	T sum = _weight[0] + _weight[1];

	return (_dx[0] * _weight[0] + _dx[1] * _weight[1]) / sum;
}

template<typename T>
void WENO3rd_Adaptive<T>::print_weights(const bool c_w, const bool s_w, const bool w_w) const {
	std::cout << "*** Upwind Left Weights _ WENO 3rd Adaptive" << std::endl << std::endl;
	print_weights(c_m, smooth_m, d_m, c_w, s_w, w_w); std::cout << std::endl;

	std::cout << "*** Upwind Right Weights _ WENO 3rd Adaptive" << std::endl << std::endl;
	print_weights(c_p, smooth_p, d_p, c_w, s_w, w_w);
}

template<typename T>
void WENO3rd_Adaptive<T>::print_weights(const T c[2][2], const T s[2][2], const T d[2], const bool c_w, const bool s_w, const bool w_w) const {
	if(c_w) {
		std::cout << "Combining Weights for 1st Derivatives" << endl;
		for(int i=1;i>=0;--i) { T sum = 0; for(size_t j=0;j<2;++j) { sum += c[i][j]; std::cout << c[i][j] << " "; } std::cout << sum << std::endl; } 
		std::cout << std::endl;
	}

	if(s_w) {
		std::cout << "Smoothness Weights" << endl;
		for(int i=1;i>=0;--i) { 
			T sum = 0;
			sum = 0; for(size_t j=0;j<2;++j) { sum += s[i][j]; std::cout << s[i][j] << " "; } 
			std::cout << sum << std::endl; 						
		} 
		std::cout << std::endl; 	
	}

	if(w_w) {
		std::cout << "WENO Weights" << endl;
		T sum = 0; for(int i=1;i>=0;--i) { sum += d[i]; std::cout << d[i] << " "; } 
		std::cout << sum << std::endl; 
	}
}

template<typename T>
void WENO3rd_Adaptive<T>::pre_compute_weights(const T _h[5]) {
	for(int i=0;i<4;++i) {
		h[i] = _h[i];
		hInv[i] = T(1.0) / h[i];
	}

	calc_combining_weight_1st_derivative(c_m, &h[0], true);
	calc_combining_weight_1st_derivative(c_p, &h[1], false);

	calc_smoothness_weight(smooth_m, &h[0]);
	calc_smoothness_weight(smooth_p, &h[1]);

	calc_weno_weight(d_m, &h[0], true);
	calc_weno_weight(d_p, &h[1], false);
}

template<typename T>
void WENO3rd_Adaptive<T>::calc_combining_weight_1st_derivative(T c[2][2], const T h[3], const bool is_positive_vel) {
	T one = T(1.0);
	T b[2][2];

	if(is_positive_vel) {
		b[1][1] = (h[0] + T(2.0) * h[1]) / ((h[0] + h[1]) * h[1]);
		b[1][0] = b[1][1] - (h[0] + h[1]) / (h[0] * h[1]);
		b[0][1] = h[1] / ((h[1] + h[2]) * h[2]);
		b[0][0] = b[0][1] - (h[1] - h[2]) / (h[1] * h[2]);
	}
	else {
		b[1][1] = h[0] / ((h[0] + h[1]) * h[1]);
		b[1][0] = b[1][1] - (h[0] - h[1]) / (h[0] * h[1]);
		b[0][1] = - h[1] / ((h[1] + h[2]) * h[2]);
		b[0][0] = b[0][1] + (h[1] + h[2]) / (h[1] * h[2]);
	}
	for(int i=0;i<2;++i)
		for(int j=0;j<2;++j) 
			c[i][j] = b[i][j] * h[1-i+j];
}

template<typename T>
void WENO3rd_Adaptive<T>::calc_smoothness_weight(T b_[2][2], const T h[4]) {
	for(int i=0;i<2;++i) {
		b_[i][1] = T(2.0) / ((h[1-i] + h[2-i]) * h[2-i]);
		b_[i][0] = b_[i][1] - T(2.0) / (h[1-i] * h[2-i]);
	}
	for(int i=0;i<2;++i) 
		for(int j=0;j<2;++j) 
			b_[i][j] *= (h[1-i+j] * h[1]);
}

template<typename T>
void WENO3rd_Adaptive<T>::calc_weno_weight(T d[2], const T h[3], const bool is_positive_vel) {
	T sum = (h[0] + h[1] + h[2]);

	if(is_positive_vel)		d[0] = (h[0] + h[1]) / sum;
	else					d[0] = h[0] / sum;
	
	d[1] = T(1.0) - d[0];
}

template<typename T>
void WENO3rd_Adaptive<T>::validate_weights_uniform_weno(const double h) {
	double h_array[4], weight[2], sum;

	for(int i=0;i<4;++i) h_array[i] = h;

	WENO3rd_Adaptive<double> weno(h_array);
	weno.print_weights();

	std::cout << "======================================================" << std::endl;
	for(int k=2;k>=0;--k) {
		sum = 0;
		combining_weight_1st_derivative_uniform(weight, k, 2, h); 
		for(int i=0;i<2;++i) {
			sum += weight[i];
			std::cout << weight[i] << " "; 
		}
		std::cout << sum << std::endl;
	}
	std::cout << "======================================================" << std::endl;

	weno3rd_weight_from_eno_uniform(weight, 2, h); sum = 0; for(int i=0;i<2;++i) { sum += weight[i]; std::cout << weight[i] << " "; } std::cout << sum << std::endl;
	weno3rd_weight_from_eno_uniform(weight, 1, h); sum = 0; for(int i=0;i<2;++i) { sum += weight[i]; std::cout << weight[i] << " "; } std::cout << sum << std::endl;

	std::cout << "======================================================" << std::endl << std::endl;
}

template<typename T>
void WENO3rd_Adaptive<T>::validate_weights_except_smoothness(const double h[4]) {
	double weight[2], sum;	

	WENO3rd_Adaptive<double> weno(h);	
	weno.print_weights();

	std::cout << "======================================================" << std::endl;
	for(int k=2;k>=0;--k) {
		sum = 0;
		combining_weight_1st_derivative(weight, k, 2, &h[2-k]); 
		for(int i=0;i<2;++i) {
			sum += weight[i];
			std::cout << weight[i] << " "; 
		}
		std::cout << sum << std::endl;
	}
	std::cout << "======================================================" << std::endl;

	weno3rd_weight_from_eno(weight, 2, &h[0]); sum = 0; for(int i=0;i<2;++i) { sum += weight[i]; std::cout << weight[i] << " "; } std::cout << sum << std::endl;
	weno3rd_weight_from_eno(weight, 1, &h[1]); sum = 0; for(int i=0;i<2;++i) { sum += weight[i]; std::cout << weight[i] << " "; } std::cout << sum << std::endl;

	std::cout << "======================================================" << std::endl << std::endl;
}

template<typename T> template<typename FUNC> 
void WENO3rd_Adaptive<T>::validate_values_uniform_weno(const double h, const FUNC& func, const double a, const double b, std::size_t n_sample) {
	basis::RandomGen<double> rand;
	WENO3rd_Adaptive<double> weno3rd_adaptive(h);

	double v[5], dx_m, dx_p;
	for(std::size_t i=0;i<n_sample;++i) {
		double x = rand.getDouble() * (b - a) + a;

		for(int k=0;k<5;++k)
			v[k] = func(x + (k-2) * h);

		weno3rd_adaptive.dx(dx_m, dx_p, v);	std::cout << dx_m << " " << dx_p << std::endl;

		dx_m = weno3rd(v, h, true);
		dx_p = weno3rd(v, h, false);		std::cout << dx_m << " " << dx_p << std::endl;

		std::cout << std::endl;
	}
}

};		// arithmetic