#ifndef _ILLUMINATION_LIBRARY_H_
#define _ILLUMINATION_LIBRARY_H_

#include <cmath>
#include "../../../Principle/Arithmetic/Vector/Vector3.h"
#include "../../IlluminationSpace.h"
#include "../Illumination/Constants.h"

namespace illuminationspace {

using arithmetic::Vector3;

class ILLUMINATIONDLL IlluminationLibrary {
public : 
	static const Vector3<real> UniformSampleSphere(real u1, real u2);
	static const real UniformSpherePdf(real theta, real phi) { return 0.5f * INV_TWOPI_i; }

	static const Vector3<real> UniformSampleHemisphere(real u1, real u2);
	static const real UniformHemispherePdf(real theta, real phi) { return INV_TWOPI_i; }

	static const Vector3<real> CosineSampleHemisphere(real u1, real u2);
	static const real CosineHemispherePdf(real costheta, real phi) { return costheta * INV_PI_i; }

	static void ConcentricSampleDisk(real u1, real u2, real *dx, real *dy);	

	static void UniformSampleTriangle(real u1, real u2, real *u, real *v);
};

};		// illuminationspace

#endif
