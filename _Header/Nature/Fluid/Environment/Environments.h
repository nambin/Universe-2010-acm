#ifndef _FLUIDSPACE_ENVIRONMENTS_H_
#define _FLUIDSPACE_ENVIRONMENTS_H_

#include <string>
#include "../../../Basis/Basis/Log.h"
#include "../../../Principle/Arithmetic/GridLibrary/Grid2.h"
#include "../../../Principle/Arithmetic/GridLibrary/Grid2x3.h"
#include "../../../Principle/Arithmetic/GridLibrary/Grid3.h"
#include "../../../Principle/Arithmetic/GridLibrary/Grid3x4.h"
#include "../../../Principle/Arithmetic/GridLibrary/GridBoundaryCondition.h"
#include "../../Nature.h"

using namespace arithmetic;

namespace fluidspace {

//#define FLUID_USE_DOUBLE_PRECISION

#ifdef FLUID_USE_DOUBLE_PRECISION
	typedef double real;
#else
	typedef float real;
#endif

typedef Grid2<bool,real>		Grid2b;
typedef Grid2<char,real>		Grid2c;
typedef Grid2<int,real>			Grid2i;
typedef Grid2<float,real>		Grid2f;
typedef Grid2<double,double>	Grid2d;
typedef Grid2<real,real>		Grid2r;

typedef Grid2x3<float>			Grid2x3f;
typedef Grid2x3<double>			Grid2x3d;
typedef Grid2x3<real>			Grid2x3r;

typedef Grid3<bool,real>		Grid3b;
typedef Grid3<char,real>		Grid3c;
typedef Grid3<int,real>			Grid3i;
typedef Grid3<float,real>		Grid3f;
typedef Grid3<double,double>	Grid3d;
typedef Grid3<real,real>		Grid3r;

typedef Grid3x4<float>			Grid3x4f;
typedef Grid3x4<double>			Grid3x4d;
typedef Grid3x4<real>			Grid3x4r;

typedef GridBoundaryCondition<real,real> GridBdry;

NATUREDLL const std::string get_concat_filename(const std::string& name, const std::size_t frame);

};		// namespace fluidspace

#endif

