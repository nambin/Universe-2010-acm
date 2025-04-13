#pragma once

namespace geometryspace {

class RungeKutta {
public : 
	RungeKutta() {}
};

class RungeKutta_Euler : public RungeKutta {
public : 
	RungeKutta_Euler() {}

	template<typename DATA, typename INTEGRATOR>
	void operator()(DATA& dst, const DATA& src, const INTEGRATOR& itg) { 
		
		if(&dst == &src) {
			basis::throwError("Error in RungeKutta_RK2::operator()");
			return;
		}

		itg(dst, src); 
	}
};

class RungeKutta_RK2 : public RungeKutta {
public : 
	RungeKutta_RK2() {}

	template<typename DATA, typename INTEGRATOR>
	void operator()(DATA& dst, const DATA& src, const INTEGRATOR& itg) { 
		typedef typename DATA::value_type T;

		if(&dst == &src) {
			basis::throwError("Error in RungeKutta_RK2::operator()");
			return;
		}

		DATA tmp; tmp.resize(dst);
		itg(dst, src); 
		itg(tmp, dst); 
		add(dst, T(0.5), tmp, T(0.5), src);
	}
};

class RungeKutta_RK3 : public RungeKutta {
public : 
	RungeKutta_RK3() {}

	template<typename DATA, typename INTEGRATOR>
	void operator()(DATA& dst, const DATA& src, const INTEGRATOR& itg) { 
		typedef typename DATA::value_type T;

		if(&dst == &src) {
			basis::throwError("Error in RungeKutta_RK3::operator()");
			return;
		}

		DATA tmp; tmp.resize(dst);
		itg(dst, src); 
		itg(tmp, dst); 
		add(dst, T(0.75), src, T(0.25), tmp);
		itg(tmp, dst);
		add(dst, T(1.0/3.0), src, T(1.0 - 1.0/3.0), tmp);
	}
};

};		// geometryspace