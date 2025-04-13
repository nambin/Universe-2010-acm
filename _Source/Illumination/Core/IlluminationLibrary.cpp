#include <stdafx.h>
#include <Illumination/IlluminationHeader/Core/IlluminationLibrary.h>
#include <Basis/Basis/BasicLibrary.h>

using namespace basis;
using namespace illuminationspace;

const Vector3<real> IlluminationLibrary::UniformSampleSphere(real u1, real u2) {
	real z = 1.0f - 2.0f * u1;
	real r = sqrt(MAX(static_cast<real>(0.0f), 1.0f - z*z));
	real phi = 2.0f * PI_i * u2;
	real x = r * cos(phi);
	real y = r * sin(phi);

	return Vector3<real>(x,y,z);
}

const Vector3<real> IlluminationLibrary::UniformSampleHemisphere(real u1, real u2) {
	real z = u1;
	real r = sqrt(MAX(static_cast<real>(0.0f), 1.0f - z*z));
	real phi = 2.0f * PI_i * u2;
	real x = r * cos(phi);
	real y = r * sin(phi);

	return Vector3<real>(x,y,z);
}

const Vector3<real> IlluminationLibrary::CosineSampleHemisphere(real u1, real u2) {
	Vector3<real> ret;
	IlluminationLibrary::ConcentricSampleDisk(u1, u2, &ret[0], &ret[1]);
	ret[2] = sqrt(MAX(static_cast<real>(0.0f), static_cast<real>(1.0f) - ret[0]*ret[0] - ret[1]*ret[1]));
	return ret;
}

void IlluminationLibrary::ConcentricSampleDisk(real u1, real u2, real *dx, real *dy) {
	real r, theta;

	// Map uniform random numbers to $[-1,1]^2$
	real sx = 2 * u1 - 1;
	real sy = 2 * u2 - 1;

	// Map square to $(r,\theta)$
	// Handle degeneracy at the origin
	if (sx == 0.0f && sy == 0.0f) {
		*dx = 0.0f;
		*dy = 0.0f;
		return;
	}
	if (sx >= -sy) {
		if (sx > sy) {	// Handle first region of disk			  
			r = sx;
			if (sy > 0.0f)  theta = sy/r;
			else			theta = 8.0f + sy/r;
		}
		else {			// Handle second region of disk			  
			r = sy;
			theta = 2.0f - sx/r;
		}
	}
	else {
		if (sx <= sy) {	  // Handle third region of disk
			r = -sx;
			theta = 4.0f - sy/r;
		}
		else {			  // Handle fourth region of disk
			r = -sy;
			theta = 6.0f + sx/r;
		}
	}
	theta *= PI_i / 4.0f;
	*dx = r*cos(theta);
	*dy = r*sin(theta);
}

void IlluminationLibrary::UniformSampleTriangle(real u1, real u2, real *u, real *v)	{
	real sqr = sqrt(u1);
	*u = 1.0f - sqr;
	*v = u2 * sqr;
}