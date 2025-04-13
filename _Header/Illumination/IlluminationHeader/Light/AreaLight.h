#ifndef _ILLUMINATION_AREALIGHT_H_
#define _ILLUMINATION_AREALIGHT_H_

#include "Light.h"
#include "../Shape/Shape.h"

namespace illuminationspace {

class ILLUMINATIONDLL AreaLight : public Light {
protected : 
	class AreaLightValue;

public : 
	AreaLight() : Light() {}
	AreaLight(const Transform &l2w, const Spectrum &le, const Shape s, int ns, int ns_stratified = 0, Sampler::SAMPLING_TYPE type = Sampler::LATINHYPERCUBE)
	: Light(new AreaLightValue(l2w, le, s, ns, ns_stratified, type)) {}

	virtual const Spectrum L(const Vector3<real> &pt, const Vector3<real> &n, const Vector3<real> &wi) const {
		return dynamic_cast<const AreaLightValue*>(value.getptr())->L(pt, n, wi);
	}

protected : 
	class ILLUMINATIONDLL AreaLightValue : public Light::LightValue {
	public : 
		AreaLightValue(const Transform &l2w, const Spectrum &le, const Shape s, int ns, int ns_stratified = 0, Sampler::SAMPLING_TYPE type = Sampler::LATINHYPERCUBE);

		AreaLightValue * clone() const { return new AreaLightValue(*this); }

		const Spectrum power(const Scene *) const	{ return _L * _area * PI_i; }
		const bool isDeltaLight() const				{ return false; }

		const real pdf(const Vector3<real> &pt, const Vector3<real> &wi) const;

		const Spectrum sample_L(const Vector3<real> &pt, Vector3<real> *wi, VisibilityTester *vis) const;
		const Spectrum sample_L(const Vector3<real> &pt, real u1, real u2, Vector3<real> *wi, real *pdf, VisibilityTester *vis) const;

		void generatePhotonPath(Vector3<real> &pos, Vector3<real> &dir) const;

		virtual const Spectrum L(const Vector3<real> &pt, const Vector3<real> &n, const Vector3<real> &wi) const {
			return dot(n,wi) > 0 ? _L : Spectrum(0.0f);
		}

	public : 
		Spectrum	_L;

		Shape		shape;
		real		_area;
	};
};

};		// illuminationspace

#endif