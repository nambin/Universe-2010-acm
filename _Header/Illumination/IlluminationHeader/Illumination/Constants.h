#ifndef _ILLUMINATION_CONSTANTS_H_
#define _ILLUMINATION_CONSTANTS_H_

#include <cmath>
#include <cfloat>
#include <cassert>
#include "../../../Basis/Basis/Constants.h"
#include "Environments.h"

namespace illuminationspace {

static const real LARGE_i		= +1.0e+6f;
static const real EPSILON_i		= +1.0e-6f;

static const real PI_i			= static_cast<real>(basis::PI);
static const real TWOPI_i		= static_cast<real>(basis::TWOPI);
static const real INV_PI_i		= static_cast<real>(basis::INV_PI);					
static const real INV_TWOPI_i	= static_cast<real>(basis::INV_TWOPI);
static const real SQRT2_i		= static_cast<real>(basis::SQRT2);
static const real SQRT3_i		= static_cast<real>(basis::SQRT3);
static const real RadToDeg_i	= static_cast<real>(basis::RadToDeg);
static const real DegToRad_i	= static_cast<real>(basis::DegToRad);

};	// namespace basis

#endif