#ifndef _ILLUMINATION_ENVIRONMENTS_H_
#define _ILLUMINATION_ENVIRONMENTS_H_

#include "../../../Basis/Basis/Environments.h"

namespace illuminationspace {

//#define ILLUMINATION_USE_DOUBLE_PRECISION

#ifdef ILLUMINATION_USE_DOUBLE_PRECISION
	typedef double real;
#else
	typedef float real;
#endif

};		// namespace illuminationspace

#endif

