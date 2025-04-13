#include <stdafx.h>
#include <Illumination/IlluminationHeader/Light/DistantLight.h>
#include <Illumination/IlluminationHeader/Scene/Scene.h>

using namespace illuminationspace;

DistantLight::DistantLightValue::DistantLightValue(const Transform &l2w, const Spectrum &radiance, const Vector3<real> &dir)
: LightValue(l2w) {
	lightDir = unitVector(LightToWorld.transformVector(dir));
	L = radiance;
}

const Spectrum DistantLight::DistantLightValue::power(const Scene *scene) const	{
	Vector3<real> center;
	real radius;
	scene->WorldBound().boundingSphere(center, radius);
	return L * PI_i * radius * radius;
}

const Spectrum DistantLight::DistantLightValue::sample_L(const Vector3<real> &pt, Vector3<real> *wi, VisibilityTester *vis) const {
	*wi = lightDir;
	vis->setRay(pt, lightDir);
	return L;
}

const Spectrum DistantLight::DistantLightValue::sample_L(const Vector3<real> &pt, real u1, real u2, Vector3<real> *wi, real *pdf, VisibilityTester *vis) const {
	*pdf = 1.0f;
	return sample_L(pt, wi, vis);
}
