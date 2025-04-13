#ifndef _ILLUMINATION_POINTLIGHT_H_
#define _ILLUMINATION_POINTLIGHT_H_

#include "Light.h"

namespace illuminationspace {

class ILLUMINATIONDLL PointLight : public Light {
protected : 
	class PointLightValue;

public : 	
	PointLight(const Transform &l2w, const Spectrum &_intensity) 
	: Light(new PointLightValue(l2w, _intensity)) {}

protected :
	class ILLUMINATIONDLL PointLightValue : public Light::LightValue {
	public : 
		PointLightValue(const Transform &l2w, const Spectrum &_intensity);	

		PointLightValue * clone() const { return new PointLightValue(*this); }

		const Spectrum power(const Scene *) const	{ return intensity * 4.0f * PI_i; }
		const bool isDeltaLight() const				{ return true; }

		const real pdf(const Vector3<real> &pt, const Vector3<real> &wi) const { return 0.0f; }

		const Spectrum sample_L(const Vector3<real> &pt, Vector3<real> *wi, VisibilityTester *vis) const;
		const Spectrum sample_L(const Vector3<real> &pt, real u1, real u2, Vector3<real> *wi, real *pdf, VisibilityTester *vis) const;

		void generatePhotonPath(Vector3<real> &pos, Vector3<real> &dir) const;

	public : 
		Vector3<real>	lightPos;
		Spectrum		intensity;
	};
};

};		// illumination

#endif