#ifndef _ILLUMINATION_MACROFACET_H_
#define _ILLUMINATION_MACROFACET_H_

#include "BxDF.h"
#include "Fresnel.h"

namespace illuminationspace {

class ILLUMINATIONDLL MicrofacetDistribution {
public : 
	virtual const real D(const Vector3<real> &wh) const = 0; 
};

class ILLUMINATIONDLL Microfacet : public BxDF {
public : 
	Microfacet(const Spectrum &r, Fresnel *f, MicrofacetDistribution *d)
	: BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)), R(r), fresnel(f), distribution(d) {}

	const Spectrum f(const Vector3<real> &wo, const Vector3<real> &wi) const;
	const real G(const Vector3<real> &wo, const Vector3<real> &wi, const Vector3<real> &wh) const;

public :
	Spectrum R;
	MicrofacetDistribution *distribution;
	Fresnel *fresnel;
};

class ILLUMINATIONDLL Blinn : public MicrofacetDistribution {
public : 
	Blinn(real e) {	exponent = e; }
	const real D(const Vector3<real> &wh) const;

public : 
	real exponent;
};

};		// illuminationspace

#endif