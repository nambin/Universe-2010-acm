#include <stdafx.h>
#include <Illumination/IlluminationHeader/Reflection/Fresnel.h>

using namespace illuminationspace;

const Spectrum Fresnel::FrDiel(real cosi, real cost, const Spectrum &etai, const Spectrum &etat) {
	Spectrum Rparl =	((etat * cosi) - (etai * cost)) /
						((etat * cosi) + (etai * cost));
	Spectrum Rperp =	((etai * cosi) - (etat * cost)) /
						((etai * cosi) + (etat * cost));
	return (Rparl*Rparl + Rperp*Rperp) * 0.5f;
}

const Spectrum Fresnel::FrCond(real cosi, const Spectrum &eta, const Spectrum &k) {
	Spectrum tmp = (eta*eta + k*k) * cosi*cosi;
	Spectrum Rparl2 =	(tmp - (2.f * eta * cosi) + Spectrum(1.0f)) /
						(tmp + (2.f * eta * cosi) + Spectrum(1.0f));
	Spectrum tmp_f = eta*eta + k*k;
	Spectrum Rperp2 =	(tmp_f - (2.f * eta * cosi) + Spectrum(cosi*cosi)) /
						(tmp_f + (2.f * eta * cosi) + Spectrum(cosi*cosi));
	return (Rparl2 + Rperp2) * 0.5f;
}

const Spectrum Fresnel::FresnelApproxEta(const Spectrum &Fr) {
	Spectrum reflectance = Fr.getClamped(0.0f, 0.999f);
	return	(Spectrum(1.0f) + reflectance.sqrt()) /
			(Spectrum(1.0f) - reflectance.sqrt());
}

const Spectrum Fresnel::FresnelApproxK(const Spectrum &Fr) {
	Spectrum reflectance = Fr.getClamped(0.0f, 0.999f);
	Spectrum sqr = (reflectance / (Spectrum(1.0f) - reflectance)); 
	return 2.0f * sqr.sqrt();
}

const Spectrum FresnelConductor::evaluate(real cosi) const {
	return Fresnel::FrCond(abs(cosi), eta, k);
}

const Spectrum FresnelDielectric::evaluate(real cosi) const {
	cosi = clamp(cosi, static_cast<real>(-1.0f), static_cast<real>(1.0f));
	bool entering = (cosi > 0.0f);
	real ei = eta_i, et = eta_t;
	if (!entering)	std::swap(ei, et);

	real eta = ei / et;

	// Compute sint^2 using Snell's law
	real sint2 = eta * eta * MAX(static_cast<real>(0.0f), 1.0f - cosi*cosi);

	if (sint2 > 1.0f)	// Total internal reflection
		return Spectrum(1.0f);
	else {
		real cost = sqrt(MAX(static_cast<real>(0.0f), 1.0f - sint2));
		return Fresnel::FrDiel(abs(cosi), cost, Spectrum(ei), Spectrum(et));
	}
}

const Spectrum FresnelDielectricSchlick::evaluate(real cosi) const {
	cosi = clamp(cosi, static_cast<real>(-1.0f), static_cast<real>(1.0f));
	bool entering = (cosi > 0.0f);
	real ei = eta_i, et = eta_t;
	if (!entering)	std::swap(ei, et);

	real eta = ei / et;

	// Compute sint^2 using Snell's law
	real sint2 = eta * eta * MAX(static_cast<real>(0.0f), 1.0f - cosi*cosi);

	if (sint2 > 1.0f)	// Total internal reflection
		return Spectrum(1.0f);
	else {
		real temp3 = 1.0f - abs(cosi);
		return Spectrum(R0 + (1.0f - R0) * temp3 * temp3 * temp3 * temp3 * temp3);
	}
}