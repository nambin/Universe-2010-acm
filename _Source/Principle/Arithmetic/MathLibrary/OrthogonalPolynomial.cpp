#include <stdafx.h>
#include <Principle/Arithmetic/MathLibrary/OrthogonalPolynomial.h>

namespace arithmetic {

void coeff_legendre(const int n, double *a, double *b) {
	for(int i=0;i<n;++i) {
		a[i] = 0.0;
		b[i] = double(i*i) / double(4*i*i-1);
	}
}

void coeff_chebyshev(const int n, double *a, double *b) {
	for(int i=0;i<n;++i) {
		a[i] = 0.0;
		b[i] = 0.25;
	}
	b[0] = 0;
}

};		// arithmetic