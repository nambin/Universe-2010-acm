#ifndef _ILLUMINATION_SPOTLIGHT_H_
#define _ILLUMINATION_SPOTLIGHT_H_

#include "Light.h"

namespace illuminationspace {

// SpotLight points down the -z axis
class ILLUMINATIONDLL SpotLight : public Light {
protected :
	class SpotLightValue;

public : 
	SpotLight(const Transform &l2w, const Spectrum &_intensity, real width, real fall)
	: Light(new SpotLightValue(l2w, _intensity, width, fall)) {}

protected : 
	class ILLUMINATIONDLL SpotLightValue : public Light::LightValue {
	public : 
		SpotLightValue(const Transform &l2w, const Spectrum &_intensity, real width, real fall);

		SpotLightValue * clone() const { return new SpotLightValue(*this); }

		const Spectrum power(const Scene *) const	{ return intensity * 2.0f * PI_i * (1.0f - 0.5f * (cosFalloffStart + cosTotalWidth)); }
		const bool isDeltaLight() const				{ return true; }

		const real pdf(const Vector3<real> &pt, const Vector3<real> &wi) const { return 0.0f; }

		const Spectrum sample_L(const Vector3<real> &pt, Vector3<real> *wi, VisibilityTester *vis) const;	
		const Spectrum sample_L(const Vector3<real> &pt, real u1, real u2, Vector3<real> *wi, real *pdf, VisibilityTester *vis) const;

		void generatePhotonPath(Vector3<real> &pos, Vector3<real> &dir) const {}		// Need to implement

		const real fallOff(const Vector3<real> &w) const;

	public : 
		real			cosTotalWidth, cosFalloffStart;
		Vector3<real>	lightPos;
		Spectrum		intensity;
	};
};

};		// illuminationspace

#endif