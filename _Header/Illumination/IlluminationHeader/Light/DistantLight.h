#ifndef _ILLUMINATION_DISTANTLIGHT_H_
#define _ILLUMINATION_DISTANTLIGHT_H_

#include "../Light/Light.h"

namespace illuminationspace {

class ILLUMINATIONDLL DistantLight : public Light {
protected : 
	class DistantLightValue;

public : 
	DistantLight(const Transform &l2w, const Spectrum &radiance, const Vector3<real> &dir)
	: Light(new DistantLightValue(l2w, radiance, dir)) {}

protected : 
	class ILLUMINATIONDLL DistantLightValue : public Light::LightValue {
	public : 
		DistantLightValue(const Transform &l2w, const Spectrum &radiance, const Vector3<real> &dir);

		DistantLightValue * clone() const { return new DistantLightValue(*this); }

		const Spectrum power(const Scene *) const;
		const bool isDeltaLight() const { return true; }

		const real pdf(const Vector3<real> &pt, const Vector3<real> &wi) const { return 0.0f; }

		const Spectrum sample_L(const Vector3<real> &pt, Vector3<real> *wi, VisibilityTester *vis) const;
		const Spectrum sample_L(const Vector3<real> &pt, real u1, real u2, Vector3<real> *wi, real *pdf, VisibilityTester *vis) const;

		void generatePhotonPath(Vector3<real> &pos, Vector3<real> &dir) const {}		// need to implement

	public : 
		Vector3<real>	lightDir;
		Spectrum		L;
	};
};

};		// illuminationspace

#endif