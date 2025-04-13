#include <stdafx.h>
#include <Principle/Arithmetic/MathLibrary/GaussQuadrature.h>
#include <Principle/Arithmetic/MathLibrary/SymmEig.h>
#include <Principle/Arithmetic/MathLibrary/GammaFunction.h>
#include <Principle/Arithmetic/MathLibrary/OrthogonalPolynomial.h>

using namespace std;

namespace arithmetic {

void gaucof(const int n, double *x, double *w, const double* const a, const double* const _b, const double amu0) {	
	double * b = new double[n];
	for (int i=0;i<n;i++)
		if (i != 0) b[i] = sqrt(_b[i]);

	SymmEig<double> sym(VectorN<double>(n,a),VectorN<double>(n,b));

	for (int i=0;i<n;i++) {
		x[i] = sym.d[i];
		w[i] = amu0 * SQ(sym.z(0,i));
	}

	delete[] b;
}

void lobatto(const int n, double *x, double *w, const double* const _a, const double* const _b, const double amu0, const double x1, const double xn) {
	if (n <= 1) {
		cout << "Error in lobatto _ n must be bigger than 1 in lobatto";
		return;
	}
	double * a = new double[n];
	double * b = new double[n];
	for(int i=0;i<n;++i) { a[i] = _a[i]; b[i] = _b[i]; }

	double det,pl,pr,p1l,p1r,pm1l,pm1r;

	pl=x1-a[0];
	pr=xn-a[0];
	pm1l=1.0;
	pm1r=1.0;
	p1l=pl;
	p1r=pr;
	for (int i=1;i<n-1;i++) {
		pl=(x1-a[i])*p1l-b[i]*pm1l;
		pr=(xn-a[i])*p1r-b[i]*pm1r;
		pm1l=p1l;
		pm1r=p1r;
		p1l=pl;
		p1r=pr;
	}
	det=pl*pm1r-pr*pm1l;
	a[n-1]=(x1*pl*pm1r-xn*pr*pm1l)/det;
	b[n-1]=(xn-x1)*pl*pr/det;
	gaucof(n,x,w,a,b,amu0);

	delete[] a; delete[] b;
}

};		// arithmetic