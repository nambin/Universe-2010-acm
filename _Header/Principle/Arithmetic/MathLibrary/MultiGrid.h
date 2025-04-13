#pragma once

#include <cmath>

namespace arithmetic {

class MultiGrid {
public : 
	template<typename T, typename GMG>
	static const T v_cycle_basis(GMG& gmg, size_t lv, const T tol, const T eps = 0.005f);

	template<typename T, typename GMG>
	static const std::size_t v_cycle(GMG& gmg, const T tol, const size_t max_iter = 200, const T eps = 0.005f);

	template<typename T, typename GMG>
	static const std::size_t v_cycle_fully_converging(GMG& gmg, size_t lv, const T tol, const size_t max_iter = 200, const T eps = 0.005);
};

template<typename T, typename GMG>
const std::size_t MultiGrid::v_cycle_fully_converging(GMG& gmg, size_t lv_des, const T tolerance, const size_t max_iter, const T eps) {

	T res = 1.0e+8;
	size_t iter = 0;

	while(res > tolerance && iter < max_iter) {
		iter++;

		gmg.smooth_ith_level(lv_des, gmg.pre_smoothing);
		res = gmg.get_residual(lv_des);

		if(res < tolerance) {
			if(lv_des == 0)
				log_info << "Error for GMG : " << res << std::endl;

			return iter + 1;
		}

		const T des_res = (lv_des == gmg.max_level()) ? tolerance : MAX(res * (1-eps), tolerance);
	
		if(lv_des != gmg.max_level()) {
			gmg.restrict_residual(lv_des);
			gmg.set_zero(lv_des+1);
			v_cycle_fully_converging(gmg, lv_des+1, tolerance, max_iter, eps);
			gmg.prolongate_correction(lv_des+1);
		}

		res = 1.0e+8; 
		while( res > des_res ) {
			gmg.smooth_ith_level(lv_des, gmg.post_smoothing);
			res = gmg.get_residual(lv_des);
		}
	}

	if(lv_des == 0)
		log_info << "Error for GMG : " << res << std::endl;

	return iter;
}

template<typename T, typename GMG>
const std::size_t MultiGrid::v_cycle(GMG& gmg, const T tolerance, const size_t maxItr, const T eps) {
	size_t iter = 0;
	T res = 1.0e+8;
	while(res > tolerance && iter < maxItr) {
		res = MultiGrid::v_cycle_basis(gmg, 0, tolerance);
		iter++;
	}
	return iter;	
}

template<typename T, typename GMG>
const T MultiGrid::v_cycle_basis(GMG& gmg, size_t lv_des, const T tolerance, const T eps) {

	gmg.smooth_ith_level(lv_des, gmg.pre_smoothing);

	T res = gmg.get_residual(lv_des);
	if(res < tolerance)
		return res;

	const T des_res = (lv_des == gmg.max_level()) ? tolerance : MAX(res * (1-eps), tolerance);
	
	if(lv_des != gmg.max_level()) {
		gmg.restrict_residual(lv_des);
		gmg.set_zero(lv_des+1);
		v_cycle_basis(gmg, lv_des+1, tolerance, eps);
		gmg.prolongate_correction(lv_des+1);
	}

	res = 1.0e+8;
	while( res > des_res ) {
		gmg.smooth_ith_level(lv_des, gmg.post_smoothing);
		res = gmg.get_residual(lv_des);
	}

	return res;
}

};		// arithmetic