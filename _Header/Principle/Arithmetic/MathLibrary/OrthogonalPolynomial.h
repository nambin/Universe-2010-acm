#pragma once

#include <cmath>
#include "../../Principle.h"

namespace arithmetic {

PRINCIPLEDLL void coeff_legendre(const int n, double *a, double *b);
PRINCIPLEDLL void coeff_chebyshev(const int n, double *a, double *b);

};		// arithmetic