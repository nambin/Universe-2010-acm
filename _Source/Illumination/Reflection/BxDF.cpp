#include <stdafx.h>
#include <Illumination/IlluminationHeader/Reflection/BxDF.h>

using namespace std;
using namespace illuminationspace;

const Spectrum BxDF::sample_f(const Vector3<real> &wo, Vector3<real> *wi, real u1, real u2, real *_pdf) const {
	*wi = IlluminationLibrary::CosineSampleHemisphere(u1, u2);
	if (wo.z() < 0.0f) (*wi)[2] *= -1.0f;
	*_pdf = pdf(wo, *wi);
	return f(wo, *wi);
}

const real BxDF::pdf(const Vector3<real> &wo, const Vector3<real> &wi) const {
	return isSameHemiSphere(wo, wi) ? abs(wi.z()) * INV_PI_i : 0.0f;
}

const Spectrum BxDF::rho(const Vector3<real> &wo, int nSamples, real *samples) const {
	cout << "No Implementation in BxDF::rho(, ,)" << endl;
	return f(wo, wo);
}

const Spectrum BxDF::rho(int nSamples, real *samples) const {
	cout << "No Implementation in BxDF::rho(,)" << endl;
	return f(Vector3<real>(0.0f), Vector3<real>(0.0f));
}