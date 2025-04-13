#ifndef _BASIS_CONSTANTS_H_
#define _BASIS_CONSTANTS_H_

#include <cmath>
#include <cfloat>
#include <cassert>
#include "Environments.h"

namespace basis {

static const double INFINITY	= FLT_MAX;

static const double LARGE		= +1.0e+6;
static const double EPSILON		= +1.0e-6;

static const double PI			= 3.14159265358979323846;		// ¥ð							
static const double TWOPI		= 2.0 * PI;
static const double INV_PI		= 0.31830988618379067154;							
static const double INV_TWOPI	= 0.15915494309189533577;
static const double SQRT2		= sqrt(2.0);
static const double SQRT3		= sqrt(3.0);
static const double RadToDeg	= 180.0 / PI;
static const double DegToRad	= PI / 180.0;

};	// namespace basis

#endif