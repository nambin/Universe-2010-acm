#include <stdafx.h>
#include <Illumination/IlluminationHeader/Reflection/Microfacet.h>

using namespace illuminationspace;

const Spectrum Microfacet::f(const Vector3<real> &wo, const Vector3<real> &wi) const {
	real cosThetaO = abs(cosTheta(wo));
	real cosThetaI = abs(cosTheta(wi));
	Vector3<real> wh = unitVector(wi + wo);
	real cosThetaH = dot(wi, wh);
	Spectrum F = fresnel->evaluate(cosThetaH);

	return R * distribution->D(wh) * G(wo, wi, wh) * F / (4.0f * cosThetaI * cosThetaO);
}

const real Microfacet::G(const Vector3<real> &wo, const Vector3<real> &wi, const Vector3<real> &wh) const {
	real NdotWh = abs(cosTheta(wh));
	real NdotWo = abs(cosTheta(wo));
	real NdotWi = abs(cosTheta(wi));
	real WOdotWH = abs(dot(wo, wh));
	return MIN3(static_cast<real>(1.0f), 2.0f * NdotWh * NdotWo / WOdotWH, 2.0f * NdotWh * NdotWi / WOdotWH);
}

const real Blinn::D(const Vector3<real> &wh) const {
	real costhetah = abs(cosTheta(wh));
	return (exponent + 2) * INV_TWOPI * pow(MAX(static_cast<real>(0.0f), costhetah), exponent);
}