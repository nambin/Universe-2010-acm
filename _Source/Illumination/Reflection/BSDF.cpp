#include <stdafx.h>
#include <Illumination/IlluminationHeader/Reflection/BSDF.h>

using namespace illuminationspace;

const Spectrum BSDF::f(const Vector3<real> &woW, const Vector3<real> &wiW, BxDFType flags) const {

	Vector3<real> wi = WorldToLocal(wiW), wo = WorldToLocal(woW);

	if(dot(wiW, n) * dot(woW, n) > 0.0f)	flags = BxDFType(flags & ~BSDF_TRANSMISSION);
	else									flags = BxDFType(flags & ~BSDF_REFLECTION);

	Spectrum f(0.0f);
	for(int i=0;i<nBxDFs;++i) 
		if(bxdfs[i]->MatchesFlags(flags))
			f += bxdfs[i]->f(wo, wi);

	return f;
}

const Spectrum BSDF::sample_f(const Vector3<real> &woW, Vector3<real> *wiW, real u1, real u2, real u3, real *pdf, BxDFType flags, BxDFType *sampledType) const {

	// Choose which _BxDF_ to sample
	int matchingComps = numComponents(flags);

	if (matchingComps == 0) {
		*pdf = 0.0f;
		return Spectrum(0.0f);
	}

	int which = MIN((int)(u3 * matchingComps), matchingComps-1);

	BxDF *bxdf = NULL;
	int count = which;
	for (int i = 0; i < nBxDFs; ++i) {
		if (bxdfs[i]->MatchesFlags(flags)) {
			if (count-- == 0) {
				bxdf = bxdfs[i];
				break;
			}
		}
	}

	// Sample chosen _BxDF_
	Vector3<real> wi;
	Vector3<real> wo = WorldToLocal(woW);
	*pdf = 0.0f;

	Spectrum f = bxdf->sample_f(wo, &wi, u1, u2, pdf);
	if (*pdf == 0.0f) return Spectrum(0.0f);
	if (sampledType) *sampledType = bxdf->type;
	*wiW = LocalToWorld(wi);

	// Compute overall PDF with all matching BxDFs
	if (!(bxdf->type & BSDF_SPECULAR) && matchingComps > 1) {
		for (int i = 0; i < nBxDFs; ++i) {
			if (bxdfs[i] != bxdf &&	bxdfs[i]->MatchesFlags(flags))
				*pdf += bxdfs[i]->pdf(wo, wi);
		}
	}
	if (matchingComps > 1) *pdf /= matchingComps;

	// Compute value of BSDF for sampled direction
	if (!(bxdf->type & BSDF_SPECULAR)) {
		f = Spectrum(0.0f);
		if (dot(*wiW, n) * dot(woW, n) > 0)	flags = BxDFType(flags & ~BSDF_TRANSMISSION);	// ignore BTDFs
		else								flags = BxDFType(flags & ~BSDF_REFLECTION);		// ignore BRDFs

		for (int i = 0; i < nBxDFs; ++i)
			if (bxdfs[i]->MatchesFlags(flags))
				f += bxdfs[i]->f(wo, wi);
	}

	return f;
}

const bool BSDF::reflectPhoton(const Vector3<real> &wIn, Vector3<real> *wOut, real u1, real u2, real u3, const Spectrum &power, Spectrum *retPower, BxDFType *reflectType) {

	int nComps = numComponents();
	if(nComps == 0) return false;

	Vector3<real> dIn = WorldToLocal(wIn), dOut;
	real pdf, avgPower;	
	Spectrum f;

	real sumPower = 0.0f;

	for(int i=0;i<nComps;++i) {

		f = bxdfs[i]->sample_f(-dIn, &dOut, u1, u2, &pdf);
		avgPower = f.sum() / 3.0f * abs(dOut.z()) / pdf ;
		sumPower += avgPower;

		if(u3 < sumPower) {
			*wOut = LocalToWorld(dOut);	
			*retPower = power * f / f.sum() * 3.0f;			
			*reflectType = bxdfs[i]->type;
			return true;
		}
	}
	return false;
}

const Spectrum BSDF::rho(BxDFType flags) const {
	Spectrum ret(0.0f);
	for (int i = 0; i < nBxDFs; ++i)
		if (bxdfs[i]->MatchesFlags(flags))
			ret += bxdfs[i]->rho();
	return ret;
}

const Spectrum BSDF::rho(const Vector3<real> &wo, BxDFType flags) const {
	Spectrum ret(0.0f);
	for (int i = 0; i < nBxDFs; ++i)
		if (bxdfs[i]->MatchesFlags(flags))
			ret += bxdfs[i]->rho(wo);
	return ret;
}