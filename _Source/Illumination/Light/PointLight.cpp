#include <stdafx.h>
#include <Illumination/IlluminationHeader/Light/PointLight.h>
#include <Illumination/IlluminationHeader/Illumination/Random.h>

using namespace illuminationspace;

PointLight::PointLightValue::PointLightValue(const Transform &l2w, const Spectrum &_intensity) 
: Light::LightValue(l2w) {
	lightPos = LightToWorld.transformPoint(Vector3<real>(0,0,0));
	intensity = _intensity;
}

const Spectrum PointLight::PointLightValue::sample_L(const Vector3<real> &pt, Vector3<real> *wi, VisibilityTester *vis) const	{
	*wi = unitVector(lightPos - pt);
	vis->setSegment(pt, lightPos);
	return intensity / sqrlength(pt, lightPos);
}

const Spectrum PointLight::PointLightValue::sample_L(const Vector3<real> &pt, real u1, real u2, Vector3<real> *wi, real *pdf, VisibilityTester *vis) const {
	*pdf = 1.0f;
	return sample_L(pt, wi, vis);
}

void PointLight::PointLightValue::generatePhotonPath(Vector3<real> &pos, Vector3<real> &dir) const {
	pos = lightPos;

	do {
		dir[0] = 2.0f * random() - 1.0f;
		dir[1] = 2.0f * random() - 1.0f;
		dir[2] = 2.0f * random() - 1.0f;
	} while(sqrmag(dir) > 1.0f);	

	normalize(dir);
}
