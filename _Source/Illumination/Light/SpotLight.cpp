#include <stdafx.h>
#include <Illumination/IlluminationHeader/Light/SpotLight.h>

using namespace illuminationspace;

SpotLight::SpotLightValue::SpotLightValue(const Transform &l2w, const Spectrum &_intensity, real width, real fall) 
: Light::LightValue(l2w) {
	lightPos = LightToWorld.transformPoint(Vector3<real>(0.0f));
	intensity = _intensity;
	cosTotalWidth = cos(Rad(width));
	cosFalloffStart = cos(Rad(fall));
}	

const Spectrum SpotLight::SpotLightValue::sample_L(const Vector3<real> &pt, Vector3<real> *wi, VisibilityTester *vis) const {
	*wi = unitVector(lightPos - pt);
	vis->setSegment(pt, lightPos);
	return intensity * fallOff(-*wi) / sqrlength(pt, lightPos);
}

const Spectrum SpotLight::SpotLightValue::sample_L(const Vector3<real> &pt, real u1, real u2, Vector3<real> *wi, real *pdf, VisibilityTester *vis) const {
	*pdf = 1.0f;
	return sample_L(pt, wi, vis);
}

// SpotLight points down the -z axis
const real SpotLight::SpotLightValue::fallOff(const Vector3<real> &w) const {
	Vector3<real> wl = unitVector(WorldToLight.transformVector(w));
	real costheta = -wl.z();
	if(costheta < cosTotalWidth)	return 0.0f;
	if(costheta > cosFalloffStart)	return 1.0f;

	real delta = (costheta - cosTotalWidth) / (cosFalloffStart - cosTotalWidth);
	return delta*delta*delta*delta;
}
