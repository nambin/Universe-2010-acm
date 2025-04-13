#pragma once

#include "Interpolation_upwind.h"

namespace arithmetic {

template<typename T>
class WENO_Adaptive {
public : 
	WENO_Adaptive() {}
	WENO_Adaptive(const T h) { set(h); }
	WENO_Adaptive(const T h[6]) { set(h); }

	void set(const T h) { const T _h[6] = {h,h,h,h,h,h}; pre_compute_weights(_h); }
	void set(const T h[6]) { pre_compute_weights(h); }

	void dx(T& dx_m, T& dx_p, const T v[7]) const;	
	void print_weights(const bool c = true, const bool s = true, const bool w = true) const;

	static void validate_weights_uniform_weno(const double h = 0.5);
	static void validate_weights_except_smoothness(const double h[6]);

	template<typename FUNC> static void validate_values_uniform_weno(const double h, const FUNC&, const double a, const double b, std::size_t n_sample = 10);
	
	const bool exportBinaryFile(std::ofstream&) const;
	const bool importBinaryFile(std::ifstream&);

private :
	void pre_compute_weights(const T h[6]);

	static void calc_combining_weight_1st_derivative(T c[3][3], const T h[5], const bool is_positive_vel);
	static void calc_smoothness_weight(T b_[3][3], T bm[3][3], T b0[3][3], T bp[3][3], const T h[5]);
	static void calc_weno_weight(T d[3], const T h[5], const bool is_positive_vel);

	const T dx(const T dev[5], const T c[3][3], const T s__[3][3], const T s_m[3][3], const T s_0[3][3], const T s_p[3][3], const T d[3]) const;
	void print_weights(const T c[3][3], const T s__[3][3], const T s_m[3][3], const T s_0[3][3], const T s_p[3][3], const T d[3], const bool c_w, const bool s, const bool w) const;

private : 
	T h[6], hInv[6];

	T c_m[3][3], c_p[3][3];

	T smooth_m__[3][3];
	T smooth_m_m[3][3], smooth_m_0[3][3], smooth_m_p[3][3];

	T smooth_p__[3][3];
	T smooth_p_m[3][3], smooth_p_0[3][3], smooth_p_p[3][3];

	T d_m[3], d_p[3];
};

template<typename T> 
const bool WENO_Adaptive<T>::exportBinaryFile(std::ofstream& ofs) const {
	ofs.write((char*)(h), sizeof(T) * 6);
	ofs.write((char*)(hInv), sizeof(T) * 6);

	ofs.write((char*)(c_m), sizeof(T) * 9);
	ofs.write((char*)(c_p), sizeof(T) * 9);

	ofs.write((char*)(smooth_m__), sizeof(T) * 9);
	ofs.write((char*)(smooth_m_m), sizeof(T) * 9);
	ofs.write((char*)(smooth_m_0), sizeof(T) * 9);
	ofs.write((char*)(smooth_m_p), sizeof(T) * 9);

	ofs.write((char*)(smooth_p__), sizeof(T) * 9);
	ofs.write((char*)(smooth_p_m), sizeof(T) * 9);
	ofs.write((char*)(smooth_p_0), sizeof(T) * 9);
	ofs.write((char*)(smooth_p_p), sizeof(T) * 9);

	ofs.write((char*)(d_m), sizeof(T) * 3);
	ofs.write((char*)(d_p), sizeof(T) * 3);

	return true;
}

template<typename T> 
const bool WENO_Adaptive<T>::importBinaryFile(std::ifstream& ifs) {
	ifs.read((char*)(h), sizeof(T) * 6);
	ifs.read((char*)(hInv), sizeof(T) * 6);

	ifs.read((char*)(c_m), sizeof(T) * 9);
	ifs.read((char*)(c_p), sizeof(T) * 9);

	ifs.read((char*)(smooth_m__), sizeof(T) * 9);
	ifs.read((char*)(smooth_m_m), sizeof(T) * 9);
	ifs.read((char*)(smooth_m_0), sizeof(T) * 9);
	ifs.read((char*)(smooth_m_p), sizeof(T) * 9);

	ifs.read((char*)(smooth_p__), sizeof(T) * 9);
	ifs.read((char*)(smooth_p_m), sizeof(T) * 9);
	ifs.read((char*)(smooth_p_0), sizeof(T) * 9);
	ifs.read((char*)(smooth_p_p), sizeof(T) * 9);

	ifs.read((char*)(d_m), sizeof(T) * 3);
	ifs.read((char*)(d_p), sizeof(T) * 3);

	return true;
}

template<typename T> 
void WENO_Adaptive<T>::dx(T& dx_m, T& dx_p, const T v[7]) const {
	T dev[6];
	for(std::size_t i=0;i<6;++i) dev[i] = (v[i+1] - v[i]) * hInv[i];

	dx_m = dx(&dev[0], c_m, smooth_m__, smooth_m_m, smooth_m_0, smooth_m_p, d_m);
	dx_p = dx(&dev[1], c_p, smooth_p__, smooth_p_m, smooth_p_0, smooth_p_p, d_p);	
}

template<typename T> 
const T WENO_Adaptive<T>::dx(const T dev[5], const T c[3][3], const T s__[3][3], const T s_m[3][3], const T s_0[3][3], const T s_p[3][3], const T d[3]) const {

	static const T sixth = (1.0 / 6.0);
	static const T eps = 1.0e-8;

	T _dx[3] = {0,0,0}, _smooth[3] = {0,0,0}, _weight[3] = {0,0,0};
	T s0, s1, s2, s3;	

	for(std::size_t i=0;i<3;++i) 
		for(std::size_t j=0;j<3;++j) 
			_dx[i] += (dev[2-i+j] * c[i][j]);

	for(std::size_t i=0;i<3;++i) {
		s0 = s1 = s2 = s3 = 0;
		for(std::size_t j=0;j<3;++j) {
			s0 += (dev[2-i+j] * s__[i][j]);
			s1 += (dev[2-i+j] * s_m[i][j]);
			s2 += (dev[2-i+j] * s_0[i][j]);
			s3 += (dev[2-i+j] * s_p[i][j]);
		}
		_smooth[i] += (s0 * s0 + sixth * (s1 * s1 + s2 * s2 + s3 * s3));
	}			

	T sum = 0;
	for(std::size_t i=0;i<3;++i) {
		_weight[i] = d[i] / SQ(eps + _smooth[i]);
		sum += _weight[i];
	}	

	return (_dx[0] * _weight[0] + _dx[1] * _weight[1] + _dx[2] * _weight[2]) / sum;
}

template<typename T>
void WENO_Adaptive<T>::pre_compute_weights(const T _h[6]) {
	for(int i=0;i<6;++i) {
		h[i] = _h[i];
		hInv[i] = T(1.0) / h[i];
	}

	calc_combining_weight_1st_derivative(c_m, &h[0], true);
	calc_combining_weight_1st_derivative(c_p, &h[1], false);

	calc_smoothness_weight(smooth_m__, smooth_m_m, smooth_m_0, smooth_m_p, &h[0]);
	calc_smoothness_weight(smooth_p__, smooth_p_m, smooth_p_0, smooth_p_p, &h[1]);

	calc_weno_weight(d_m, &h[0], true);
	calc_weno_weight(d_p, &h[1], false);
}

template<typename T>
void WENO_Adaptive<T>::calc_combining_weight_1st_derivative(T c[3][3], const T h[5], const bool is_positive_vel) {	
	T one = T(1.0);
	T b[3][3];

	if(is_positive_vel) {
		b[2][2] = one / (h[0] + h[1] + h[2]) + one / (h[1] + h[2]) + one / h[2];
		b[2][1] = b[2][2] - (h[0] + h[1] + h[2]) * (h[1] + h[2]) / ((h[0] + h[1]) * h[1] * h[2]);
		b[2][0] = b[2][1] + (h[0] + h[1] + h[2]) * h[2] / (h[0] * h[1] * (h[1] + h[2]));

		b[1][2] = (h[1] + h[2]) * h[2] / ((h[1] + h[2] + h[3]) * (h[2] + h[3]) * h[3]);
		b[1][1] = b[1][2] + one / (h[1] + h[2]) + one / h[2] - one / h[3];
		b[1][0] = b[1][1] - ((h[1] + h[2]) * h[3] / (h[1] * h[2] * (h[2] + h[3])));

		b[0][2] = - h[2] * h[3] / ((h[2] + h[3] + h[4]) * (h[3] + h[4]) * h[4]);
		b[0][1] = b[0][2] + (h[2] * (h[3] + h[4]) / ((h[2] + h[3]) * h[3] * h[4]));
		b[0][0] = b[0][1] + one / h[2] - one / h[3] - one / (h[3] + h[4]);
	}
	else {
		b[2][2] = (h[0] + h[1]) * h[1] / ((h[0] + h[1] + h[2]) * (h[1] + h[2]) * h[2]);
		b[2][1] = b[2][2] + one / (h[0] + h[1]) + one / h[1] - one / h[2];
		b[2][0] = b[2][1] - ((h[0] + h[1]) * h[2] / (h[0] * h[1] * (h[1] + h[2])));

		b[1][2] = - h[1] * h[2] / ((h[1] + h[2] + h[3]) * (h[2] + h[3]) * h[3]);
		b[1][1] = b[1][2] + (h[1] * (h[2] + h[3]) / ((h[1] + h[2]) * h[2] * h[3]));
		b[1][0] = b[1][1] + one / h[1] - one / h[2] - one / (h[2] + h[3]);

		b[0][2] = h[2] * (h[2] + h[3]) / ((h[2] + h[3] + h[4]) * (h[3] + h[4]) * h[4]);
		b[0][1] = b[0][2] - h[2] * (h[2] + h[3] + h[4]) / ((h[2] + h[3]) * h[3] * h[4]);
		b[0][0] = b[0][1] + (h[2] + h[3]) * (h[2] + h[3] + h[4]) / (h[2] * h[3] * (h[3] + h[4]));
	}
	for(int i=0;i<3;++i)
		for(int j=0;j<3;++j) 
			c[i][j] = b[i][j] * h[2-i+j];
}

template<typename T>
void WENO_Adaptive<T>::calc_smoothness_weight(T b__[3][3], T bm[3][3], T b0[3][3], T bp[3][3], const T h[5]) {
	for(int i=0;i<3;++i) {
		b__[i][2] = T(6.0) / ((h[2-i] + h[3-i] + h[4-i]) * (h[3-i] + h[4-i]) * h[4-i]);
		b__[i][1] = b__[i][2] - T(6.0) / ((h[2-i] + h[3-i]) * h[3-i] * h[4-i]);
		b__[i][0] = b__[i][1] + T(6.0) / (h[2-i] * h[3-i] * (h[3-i] + h[4-i]));
	}

	T two = T(2.0), four = T(4.0), half = T(0.5);

	bm[2][2] = two * (h[0] + two * h[1]) / ((h[0] + h[1] + h[2]) * (h[1] + h[2]) * h[2]);
	bm[2][1] = bm[2][2] - two * (h[0] + two * h[1] - h[2]) / ((h[0] + h[1]) * h[1] * h[2]);
	bm[2][0] = bm[2][1] + two * (h[0] + h[1] - h[2]) / (h[0] * h[1] * (h[1] + h[2]));

	bm[1][2] = two * (h[1] - h[2]) / ((h[1] + h[2] + h[3]) * (h[2] + h[3]) * h[3]);
	bm[1][1] = bm[1][2] - two * (h[1] - h[2] - h[3]) / ((h[1] + h[2]) * h[2] * h[3]);
	bm[1][0] = bm[1][1] + two * (h[1] - two * h[2] - h[3]) / (h[1] * h[2] * (h[2] + h[3]));

	bm[0][2] = - (four * h[2] + two * h[3]) / ((h[2] + h[3] + h[4]) * (h[3] + h[4]) * h[4]);
	bm[0][1] = bm[0][2] + two * (two * h[2] + h[3] + h[4]) / ((h[2] + h[3]) * h[3] * h[4]);
	bm[0][0] = bm[0][1] - two * (two * h[2] + two * h[3] + h[4]) / (h[2] * h[3] * (h[3] + h[4]));

	for(int i=0;i<3;++i) {
		for(int j=0;j<3;++j) {
			b0[i][j] = bm[i][j] + half * h[2] * b__[i][j];
			bp[i][j] = bm[i][j] + h[2] * b__[i][j];
		}		
	}

	for(int i=0;i<3;++i) {
		for(int j=0;j<3;++j) {
			b__[i][j] *= (h[2-i+j] * h[2] * h[2]);

			bm[i][j] *= (T(1.0) * h[2-i+j] * h[2]);
			b0[i][j] *= (T(2.0) * h[2-i+j] * h[2]);
			bp[i][j] *= (T(1.0) * h[2-i+j] * h[2]);
		}
	}
}

template<typename T>
void WENO_Adaptive<T>::calc_weno_weight(T d[3], const T h[5], const bool is_positive_vel) {
	T sum = (h[0] + h[1] + h[2] + h[3] + h[4]);

	if(is_positive_vel) {
		d[2] = h[3] * (h[3] + h[4]) / ((h[0] + h[1] + h[2] + h[3]) * sum);
		d[1] = (h[0] + h[1] + h[2]) * (h[3] + h[4]) * (sum + h[1] + h[2] + h[3]) / ((sum - h[4]) * (sum - h[0]) * sum);
		d[0] = (h[1] + h[2]) * (h[0] + h[1] + h[2]) / ((sum - h[0]) * sum);
	}
	else {
		d[2] = (h[2] + h[3]) * (h[2] + h[3] + h[4]) / ((sum - h[4]) * sum);
		d[1] = (h[0] + h[1]) * (h[2] + h[3] + h[4]) * (sum + h[1] + h[2] + h[3]) / ((sum - h[4]) * (sum - h[0]) * sum);
		d[0] = h[1] * (h[0] + h[1]) / ((sum - h[0]) * sum);
	}
}

template<typename T>
void WENO_Adaptive<T>::print_weights(const T c[3][3], const T s__[3][3], const T s_m[3][3], const T s_0[3][3], const T s_p[3][3], const T d[3], const bool c_w, const bool s_w, const bool w_w) const {
	if(c_w) {
		std::cout << "Combining Weights for 1st Derivatives" << endl;
		for(int i=2;i>=0;--i) { T sum = 0; for(size_t j=0;j<3;++j) { sum += c[i][j]; std::cout << c[i][j] << " "; } std::cout << sum << std::endl; } std::cout << std::endl;
	}

	if(s_w) {
		std::cout << "Smoothness Weights" << endl;
		for(int i=2;i>=0;--i) { 
			T sum = 0;
			sum = 0; for(size_t j=0;j<3;++j) { sum += s__[i][j]; std::cout << s__[i][j] << " "; } std::cout << sum << std::endl; 
			sum = 0; for(size_t j=0;j<3;++j) { sum += s_m[i][j]; std::cout << s_m[i][j] << " "; } std::cout << sum << std::endl; 
			sum = 0; for(size_t j=0;j<3;++j) { sum += s_0[i][j]; std::cout << s_0[i][j] << " "; } std::cout << sum << std::endl; 
			sum = 0; for(size_t j=0;j<3;++j) { sum += s_p[i][j]; std::cout << s_p[i][j] << " "; } std::cout << sum << std::endl; 
			std::cout << std::endl;
		} 
	}

	if(w_w) {
		std::cout << "WENO Weights" << endl;
		T sum = 0; for(int i=2;i>=0;--i) { sum += d[i]; std::cout << d[i] << " "; } std::cout << sum << std::endl; 
	}
}

template<typename T>
void WENO_Adaptive<T>::print_weights(const bool c_w, const bool s_w, const bool w_w) const {
	std::cout << "*** Upwind Left Weights _ WENO Adaptive" << std::endl << std::endl;
	print_weights(c_m, smooth_m__, smooth_m_m, smooth_m_0, smooth_m_p, d_m, c_w, s_w, w_w); std::cout << std::endl;

	std::cout << "*** Upwind Right Weights _ WENO Adaptive" << std::endl << std::endl;
	print_weights(c_p, smooth_p__, smooth_p_m, smooth_p_0, smooth_p_p, d_p, c_w, s_w, w_w);
}

template<typename T>
void WENO_Adaptive<T>::validate_weights_uniform_weno(const double h) {	
	double h_array[6], weight[3], sum;

	for(int i=0;i<6;++i) h_array[i] = h;

	WENO_Adaptive<double> weno(h_array);
	weno.print_weights();

	std::cout << "======================================================" << std::endl;
	for(int k=3;k>=0;--k) {
		sum = 0;
		combining_weight_1st_derivative_uniform(weight, k, 3, h); 
		for(int i=0;i<3;++i) {
			sum += weight[i];
			std::cout << weight[i] << " "; 
		}
		std::cout << sum << std::endl;
	}
	std::cout << "======================================================" << std::endl;

	weno_weight_from_eno_uniform(weight, 3, h); sum = 0; for(int i=0;i<3;++i) { sum += weight[i]; std::cout << weight[i] << " "; } std::cout << sum << std::endl;
	weno_weight_from_eno_uniform(weight, 2, h); sum = 0; for(int i=0;i<3;++i) { sum += weight[i]; std::cout << weight[i] << " "; } std::cout << sum << std::endl;

	std::cout << "======================================================" << std::endl << std::endl;
}

template<typename T>
void WENO_Adaptive<T>::validate_weights_except_smoothness(const double h[6]) {	
	double weight[3], sum;	

	WENO_Adaptive<double> weno(h);
	//weno.print_weights(true, false, true);
	weno.print_weights();

	std::cout << "======================================================" << std::endl;
	for(int k=3;k>=0;--k) {
		sum = 0;
		combining_weight_1st_derivative(weight, k, 3, &h[3-k]); 
		for(int i=0;i<3;++i) {
			sum += weight[i];
			std::cout << weight[i] << " "; 
		}
		std::cout << sum << std::endl;
	}
	std::cout << "======================================================" << std::endl;

	weno_weight_from_eno(weight, 3, &h[0]); sum = 0; for(int i=0;i<3;++i) { sum += weight[i]; std::cout << weight[i] << " "; } std::cout << sum << std::endl;
	weno_weight_from_eno(weight, 2, &h[1]); sum = 0; for(int i=0;i<3;++i) { sum += weight[i]; std::cout << weight[i] << " "; } std::cout << sum << std::endl;

	std::cout << "======================================================" << std::endl << std::endl;
}

template<typename T> template<typename FUNC> 
void WENO_Adaptive<T>::validate_values_uniform_weno(const double h, const FUNC& func, const double a, const double b, std::size_t n_sample) {
	basis::RandomGen<double> rand;
	WENO_Adaptive<double> weno_adaptive(h);

	double v[7], dx_m, dx_p;
	for(std::size_t i=0;i<n_sample;++i) {
		double x = rand.getDouble() * (b - a) + a;

		for(int k=0;k<7;++k)
			v[k] = func(x + (k-3) * h);

		weno_adaptive.dx(dx_m, dx_p, v);	std::cout << dx_m << " " << dx_p << std::endl;

		dx_m = weno(v, h, true);
		dx_p = weno(v, h, false);			std::cout << dx_m << " " << dx_p << std::endl;

		std::cout << std::endl;
	}
}

};		// arithmetic