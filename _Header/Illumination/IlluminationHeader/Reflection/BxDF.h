#ifndef _ILLUMINATION_BXDF_H_
#define _ILLUMINATION_BXDF_H_

#include "../Core/IlluminationLibrary.h"
#include "../Core/Spectrum.h"

namespace illuminationspace {

inline real cosTheta(const Vector3<real> &v)	{ return v.z(); }
inline real sinTheta(const Vector3<real> &v)	{ return sqrt(MAX(static_cast<real>(0.0f), 1.0f-v.z()*v.z())); }
inline real sinTheta2(const Vector3<real> &v)	{ return 1.0f - v.z()*v.z(); }
inline real cosPhi(const Vector3<real> &v)		{ return v.x() / sinTheta(v); }
inline real sinPhi(const Vector3<real> &v)		{ return v.y() / sinTheta(v); }
inline bool isSameHemiSphere(const Vector3<real> &w, const Vector3<real> &wp) { return w.z() * wp.z() > 0.0f; }

enum BxDFType {
	BSDF_REFLECTION			= 1<<0,
	BSDF_TRANSMISSION		= 1<<1,
	BSDF_DIFFUSE			= 1<<2,
	BSDF_GLOSSY				= 1<<3,
	BSDF_SPECULAR			= 1<<4,
	BSDF_ALL_TYPES			= BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR,
	BSDF_ALL_REFLECTION		= BSDF_REFLECTION | BSDF_ALL_TYPES,
	BSDF_ALL_TRANSMISSION	= BSDF_TRANSMISSION | BSDF_ALL_TYPES,
	BSDF_ALL				= BSDF_ALL_REFLECTION | BSDF_ALL_TRANSMISSION
};

class ILLUMINATIONDLL BxDF {
public : 
	BxDF(BxDFType t) : type(t) {}

	const bool MatchesFlags(BxDFType flags) const {	return (type & flags) == type;	}

	virtual const Spectrum f(const Vector3<real> &wo, const Vector3<real> &wi) const = 0;
	virtual const Spectrum sample_f(const Vector3<real> &wo, Vector3<real> *wi, real u1, real u2, real *_pdf) const;	// sample direction of reflection

	virtual const real pdf(const Vector3<real> &wo, const Vector3<real> &wi) const;

	virtual const Spectrum rho(const Vector3<real> &wo, int nSamples = 16, real *samples = NULL) const;		// Total Power of reflection
	virtual const Spectrum rho(int nSamples = 16, real *samples = NULL) const;								// Total Power of reflection

public : 
	const BxDFType type;
};

};		// illuminationspace

#endif