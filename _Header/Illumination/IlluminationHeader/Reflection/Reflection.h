#ifndef _ILLUMINATION_REFLECTION_H_
#define _ILLUMINATION_REFLECTION_H_

#include "BxDF.h"
#include "Fresnel.h"
//#include "Microfacet.h"

namespace illuminationspace {

class SpecularReflection : public BxDF {
public : 
	SpecularReflection(const Spectrum &r, Fresnel *f) 
	: BxDF(BxDFType(BSDF_REFLECTION | BSDF_SPECULAR)), R(r), fresnel(f) {}

	const Spectrum f(const Vector3<real> &wo, const Vector3<real> &wi) const {
		return Spectrum(0.0f);
	}

	const Spectrum sample_f(const Vector3<real> &wo, Vector3<real> *wi, real u1, real u2, real *pdf) const;

	const real pdf(const Vector3<real> &wo, const Vector3<real> &wi) const {
		return 0.0f;
	}

public :
	Spectrum	R;
	Fresnel*	fresnel;
};

class SpecularTransmission : public BxDF {
public :
	SpecularTransmission(const Spectrum &t, real ei, real et, FresnelDielectric *f)
	: BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR)), T(t), eta_i(ei), eta_t(et), fresnel(f) {}

	const Spectrum f(const Vector3<real> &wo, const Vector3<real> &wi) const {
		return Spectrum(0.0f);
	}

	const Spectrum sample_f(const Vector3<real> &wo, Vector3<real> *wi, real u1, real u2, real *pdf) const;

	const real pdf(const Vector3<real> &wo, const Vector3<real> &wi) const {
		return 0.0f;
	}	

public :
	Spectrum				T;
	real					eta_i, eta_t;
	FresnelDielectric*		fresnel;
};

class Lambertian : public BxDF {
public :
	Lambertian(const Spectrum& reflectance)
		: BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R(reflectance), R_PI(reflectance * INV_PI_i) {};

	const Spectrum f(const Vector3<real> &wo, const Vector3<real> &wi) const {
		return R_PI;
	}

	const Spectrum rho(int, real *)							{ return R; }
	const Spectrum rho(const Vector3<real>&, int, real *)	{ return R; }

public :
	Spectrum R, R_PI;
};

class OrenNayar : public BxDF {
public : 
	OrenNayar(const Spectrum &r, real sigma) : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R(r), R_PI(r * INV_PI_i) {
		real sig = Rad(sigma);
		real sig2 = sig * sig;
		A = 1.0f - (sig2 / (2.0f * (sig2 + 0.33f)));
		B = 0.45f * sig2 / (sig2 + 0.09f);
	}
	const Spectrum f(const Vector3<real> &wo, const Vector3<real> &wi) const;

public :
	Spectrum	R, R_PI;
	real		A,B;
};

//class FresnelBlend : public BxDF {
//public :
//	FresnelBlend(const Spectrum &d, const Spectrum &s, MicrofacetDistribution *dist)
//		: BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)), Rd(d), Rs(s), distribution(dist) {};
//
//	Spectrum SchlickFresnel(real costheta) const {
//		real tmp = 1.0f - costheta;
//		return Rs + (Spectrum(1.0f) - Rs) * tmp * tmp * tmp * tmp * tmp;
//	}
//	Spectrum f(const Vector3<real> &wo, const Vector3<real> &wi) const;
//
//public :
//	Spectrum Rd, Rs;
//	MicrofacetDistribution *distribution;
//};

};		// illuminationspace

#endif
