#ifndef _ARITHMETIC_CUBIC_EQUATION_SOLVER_CARDANO_H_
#define _ARITHMETIC_CUBIC_EQUATION_SOLVER_CARDANO_H_

#include "Complex.h"

namespace arithmetic {

template<typename T>
class CubicEquationSolver_Cardano {
public : 
	CubicEquationSolver_Cardano();
	void solve(T a, T b, T c, T d);		// a * x^3 + b * x^2 + c * x + d = 0

	template <typename T> friend std::ostream& operator<<(std::ostream&, const CubicEquationSolver_Cardano<T>&);

private : 
	const T cubeRoot(T x) const;

	T inv2, inv3, inv4, inv27, rootThree;

public : 
	short		numReal;
	Complex<T>	sol[3];
};

template<typename T>
CubicEquationSolver_Cardano<T>::CubicEquationSolver_Cardano() {
	inv2 = T(1.0 / 2.0); 
	inv3 = T(1.0 / 3.0); 
	inv4 = T(1.0 / 4.0); 
	inv27 = T(1.0 / 27.0); 
	rootThree = T(sqrt(3.0));
}

template<typename T>
const T CubicEquationSolver_Cardano<T>::cubeRoot(T x) const {
	if (x < 0)	return - ::pow(-x, inv3);
	else		return ::pow(x, inv3);
}

template<typename T>
void CubicEquationSolver_Cardano<T>::solve(T a, T b, T c, T d) {
	
	// find the discriminant
	T a2 = a*a;
	T a3 = a2*a;
	T b2 = b*b; 
	T b3 = b2*b;	
	
	T f = (3 * c / a - b2 / a2) * inv3;
	T g = (2 * b3 / a3 - 9 * b * c / a2 + 27 * d / a) * inv27;

	T g2 = g*g, f3 = f*f*f;

	T h = g2 * inv4 + f3 * inv27;
		
	if (f == 0 && g == 0 && h == 0) {	// 3 equal roots
										// when f, g, and h all equal 0 the roots can be found by the following line
		sol[0] = sol[1] = sol[2] = - cubeRoot(d / a);		
		numReal = 1;
	}
	else if (h <= 0) {		// 3 real roots
							// complicated maths making use of the method
		
		T i = ::pow(g2 * inv4 - h, inv2);
		T j = cubeRoot(i);
		T k = acos(-(g * inv2 / i));
		T m = cos(k * inv3);
		T n = rootThree * sin(k * inv3);
		T p = - (b * inv3 / a);

		sol[0].set(2 * j * m + p, 0);
		sol[1].set(-j * (m + n) + p, 0);
		sol[2].set(-j * (m - n) + p, 0);
		numReal = 3;
	}
	else if (h > 0) {		// 1 real root and 2 complex roots
							// complicated maths making use of the method
		
		T r = -(g * inv2) + ::pow(h, inv2);
		T s = cubeRoot(r);
		T t = -(g * inv2) - ::pow(h, inv2);
		T u = cubeRoot(t);
		T p = -(b / (3 * a));

		sol[0].set((s + u) + p, 0);
		sol[1].set(-(s + u) * inv2 + p, (s - u) * rootThree * inv2);
		sol[2].set(-(s + u) * inv2 + p, - (s - u) * rootThree * inv2);
		numReal = 1;
	}
}

template<typename T> std::ostream& operator<<(std::ostream& os, const CubicEquationSolver_Cardano<T>& solver) {
	os << solver.sol[0] << std::endl;
	os << solver.sol[1] << std::endl;
	os << solver.sol[2] << std::endl;
	return os;
}

};		// arithmetic

#endif