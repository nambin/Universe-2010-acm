#include <stdafx.h>
#include <Principle/Arithmetic/MathLibrary/GaussQuadrature.h>
#include <Principle/Arithmetic/MathLibrary/GammaFunction.h>
#include <Principle/Arithmetic/MathLibrary/OrthogonalPolynomial.h>

using namespace std;

namespace arithmetic {

void gauss_reverse_order(const int n, double *x, double *w) {
	for(int k=0;k<(n/2);++k) {
		double tmp = x[k]; x[k] = x[n-1-k]; x[n-1-k] = tmp;
		tmp = w[k]; w[k] = w[n-1-k]; w[n-1-k] = tmp;
	}
}

void gauss_normalize(const int n, double *x, double *w, const double x1, const double x2) {
	double diff = x2 - x1; 
	double _w = diff / 2.0;

	for(int k=0;k<n;++k) {
		x[k] = (x[k] + 1) * 0.5 * diff + x1;
		w[k] *= _w;
	}
}

void gauss_legendre(const int n, double *x, double *w, const double x1, const double x2) {
	const double EPS = 1.0e-14;
	double z1,z,xm,xl,pp,p3,p2,p1;

	int m=(n+1)/2;

	xm = 0.5*(x2+x1);
	xl = 0.5*(x2-x1);

	for (int i=0;i<m;i++) {
		z = cos(3.141592654 * (i+0.75) / (n+0.5));
		do {
			p1=1.0;
			p2=0.0;
			for (int j=0;j<n;j++) {
				p3=p2;
				p2=p1;
				p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
			}
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
		} while (::abs(z-z1) > EPS);

		x[i]=xm-xl*z;
		x[n-1-i]=xm+xl*z;
		w[i]=2.0*xl/((1.0-z*z)*pp*pp);
		w[n-1-i]=w[i];
	}
}

void gauss_chebyshev(const int n, double *x, double *w, const double x1, const double x2) {
	//gauss_jacobi(n, x, w, -0.5, -0.5, x1, x2);
	for(int i=0;i<n;++i) {
		x[i] = cos(double(PI * (i + 0.5)) / n);
		w[i] = PI / double(n);
	}

	gauss_reverse_order(n, x, w);
	gauss_normalize(n, x, w, x1, x2);
}

void gauss_jacobi(const int n, double *x, double *w, const double alf, const double bet, const double x1, const double x2) {
	const int MAXIT=10;
	const double EPS=1.0e-14;
	int i,its,j;
	double alfbet,an,bn,r1,r2,r3;
	double a,b,c,p1,p2,p3,pp,temp,z,z1;
	
	for (i=0;i<n;i++) {
		if (i == 0) {
			an=alf/n;
			bn=bet/n;
			r1=(1.0+alf)*(2.78/(4.0+n*n)+0.768*an/n);
			r2=1.0+1.48*an+0.96*bn+0.452*an*an+0.83*an*bn;
			z=1.0-r1/r2;
		} else if (i == 1) {
			r1=(4.1+alf)/((1.0+alf)*(1.0+0.156*alf));
			r2=1.0+0.06*(n-8.0)*(1.0+0.12*alf)/n;
			r3=1.0+0.012*bet*(1.0+0.25*::abs(alf))/n;
			z -= (1.0-z)*r1*r2*r3;
		} else if (i == 2) {
			r1=(1.67+0.28*alf)/(1.0+0.37*alf);
			r2=1.0+0.22*(n-8.0)/n;
			r3=1.0+8.0*bet/((6.28+bet)*n*n);
			z -= (x[0]-z)*r1*r2*r3;
		} else if (i == n-2) {
			r1=(1.0+0.235*bet)/(0.766+0.119*bet);
			r2=1.0/(1.0+0.639*(n-4.0)/(1.0+0.71*(n-4.0)));
			r3=1.0/(1.0+20.0*alf/((7.5+alf)*n*n));
			z += (z-x[n-4])*r1*r2*r3;
		} else if (i == n-1) {
			r1=(1.0+0.37*bet)/(1.67+0.28*bet);
			r2=1.0/(1.0+0.22*(n-8.0)/n);
			r3=1.0/(1.0+8.0*alf/((6.28+alf)*n*n));
			z += (z-x[n-3])*r1*r2*r3;
		} else {
			z=3.0*x[i-1]-3.0*x[i-2]+x[i-3];
		}
		alfbet=alf+bet;
		for (its=1;its<=MAXIT;its++) {
			temp=2.0+alfbet;
			p1=(alf-bet+temp*z)/2.0;
			p2=1.0;
			for (j=2;j<=n;j++) {
				p3=p2;
				p2=p1;
				temp=2*j+alfbet;
				a=2*j*(j+alfbet)*(temp-2.0);
				b=(temp-1.0)*(alf*alf-bet*bet+temp*(temp-2.0)*z);
				c=2.0*(j-1+alf)*(j-1+bet)*temp;
				p1=(b*p2-c*p3)/a;
			}
			pp=(n*(alf-bet-temp*z)*p1+2.0*(n+alf)*(n+bet)*p2)/(temp*(1.0-z*z));
			z1=z;
			z=z1-p1/pp;
			if (::abs(z-z1) <= EPS) break;
		}
		if (its > MAXIT) {
			cout << ("too many iterations in gaujac") << endl;
			exit(0);
		}
		x[i]=z;
		w[i]=::exp(gammln(alf+n)+gammln(bet+n)-gammln(n+1.0)-gammln(n+alfbet+1.0))*temp*::pow(2.0,alfbet)/(pp*p2);
	}

	// make reverse order
	gauss_reverse_order(n, x, w);
	gauss_normalize(n, x, w, x1, x2);
}

void gauss_legendre_using_gaucof(const int n, double *x, double *w, const double x1, const double x2) {
	double * a = new double[n];
	double * b = new double[n];

	coeff_legendre(n, a, b);
	gaucof(n, x, w, a, b, 2.0);

	gauss_reverse_order(n, x, w);
	gauss_normalize(n, x, w, x1, x2);

	delete[] a; delete[] b;
}

// Wrong Result, don't know why
void gauss_chebyshev_using_gaucof(const int n, double *x, double *w, const double x1, const double x2) {
	double * a = new double[n];
	double * b = new double[n];

	coeff_chebyshev(n, a, b);
	gaucof(n, x, w, a, b, PI);

	gauss_reverse_order(n, x, w);
	gauss_normalize(n, x, w, x1, x2);

	delete[] a; delete[] b;
}

void gauss_lobatto_legendre(const int n, double *x, double *w, const double x1, const double x2) {
	double * a = new double[n];
	double * b = new double[n];

	coeff_legendre(n, a, b);
	lobatto(n, x, w, a, b, 2.0, -1, 1);

	gauss_reverse_order(n, x, w);
	gauss_normalize(n, x, w, x1, x2);

	delete[] a; delete[] b;
}

};		// arithmetic