#pragma once

#include <cmath>
#include "../../Principle.h"

namespace arithmetic {

PRINCIPLEDLL void gauss_legendre(const int n, double *x, double *w, const double x1 = -1.0, const double x2 = 1.0);
PRINCIPLEDLL void gauss_chebyshev(const int n, double *x, double *w, const double x1 = -1.0, const double x2 = 1.0);
PRINCIPLEDLL void gauss_jacobi(const int n, double *x, double *w, const double alf, const double bet, const double x1 = -1.0, const double x2 = 1.0);

PRINCIPLEDLL void gauss_legendre_using_gaucof(const int n, double *x, double *w, const double x1 = -1.0, const double x2 = 1.0);
PRINCIPLEDLL void gauss_chebyshev_using_gaucof(const int n, double *x, double *w, const double x1 = -1.0, const double x2 = 1.0);	// Wrong Result, don't know why

PRINCIPLEDLL void gauss_lobatto_legendre(const int n, double *x, double *w, const double x1 = -1.0, const double x2 = 1.0);

PRINCIPLEDLL void gaucof(const int n, double *x, double *w, const double* const a, const double* const b, const double amu0);
PRINCIPLEDLL void lobatto(const int n, double *x, double *w, const double* const a, const double* const b, const double amu0, const double x1, const double xn);


};		// arithmetic