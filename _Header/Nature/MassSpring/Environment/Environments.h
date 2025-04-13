#ifndef _MASS_SPRING_SPACE_ENVIRONMENTS_H_
#define _MASS_SPRING_SPACE_ENVIRONMENTS_H_

namespace springspace {

#define SPRING_USE_DOUBLE_PRECISION

#ifdef SPRING_USE_DOUBLE_PRECISION
	typedef double real;
#else
	typedef float real;
#endif

};		// springspace

#endif