#include <stdafx.h>
#include <Illumination/IlluminationHeader/Reflection/Reflection.h>

using namespace illuminationspace;

const Spectrum SpecularReflection::sample_f(const Vector3<real> &wo, Vector3<real> *wi, real u1, real u2, real *pdf) const {
	*wi = Vector3<real>(-wo.x(), -wo.y(), wo.z());
	*pdf = 1.0f;
	return R * fresnel->evaluate(cosTheta(wo)) / abs(cosTheta(*wi));		
}

const Spectrum SpecularTransmission::sample_f(const Vector3<real> &wo, Vector3<real> *wi, real u1, real u2, real *pdf) const {
	// Figure out which eta is incident and which is transmitted
	bool entering = (cosTheta(wo) > 0.0f);
	real ei = eta_i, et = eta_t;
	if (!entering)
		std::swap(ei, et);

	// Compute transmitted ray direction
	real sini2 = sinTheta2(wo);
	real eta = ei / et;
	real sint2 = eta * eta * sini2;

	if (sint2 > 1.0f) return Spectrum(0.0f); // Total internal reflection for transmission

	real cost = sqrt(MAX(static_cast<real>(0.0f), 1.0f - sint2));
	if (entering) cost = -cost;
	real sintOverSini = eta;

	*wi = Vector3<real>(sintOverSini * -wo.x(), sintOverSini * -wo.y(), cost);
	*pdf = 1.0f;

	return T * (Spectrum(1.0f) - fresnel->evaluate(cosTheta(wo))) / abs(cosTheta(*wi));		
}

const Spectrum OrenNayar::f(const Vector3<real> &wo, const Vector3<real> &wi) const {
	real sinthetai = sinTheta(wi);
	real sinthetao = sinTheta(wo);

	// Compute cosine term of OrenNayar model
	real sinphii = sinPhi(wi), cosphii = cosPhi(wi);
	real sinphio = sinPhi(wo), cosphio = cosPhi(wo);
	real dcos = cosphii * cosphio + sinphii * sinphio;

	real maxcos = MAX(static_cast<real>(0.0f), dcos);

	// Compute sine and tangent terms of OrenNayar model
	real sinalpha, tanbeta;
	if (abs(cosTheta(wi)) > abs(cosTheta(wo))) {
		sinalpha = sinthetao;
		tanbeta = sinthetai / abs(cosTheta(wi));
	}
	else {
		sinalpha = sinthetai;
		tanbeta = sinthetao / abs(cosTheta(wo));
	}
	return R_PI * (A + B * maxcos * sinalpha * tanbeta);
}