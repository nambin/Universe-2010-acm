#include <stdafx.h>
#include <Illumination/IlluminationHeader/Light/Light.h>
#include <Illumination/IlluminationHeader/Scene/Scene.h>

using namespace illuminationspace;

Light::LightValue::LightValue() 
: LightToWorld(NULL), WorldToLight(NULL), nSamples(1), nSamples_stratified(0), samplingType(Sampler::RANDOM) {
}

Light::LightValue::LightValue(const Transform &l2w, int ns, int ns_stratified, Sampler::SAMPLING_TYPE type) 
: nSamples(MAX(1,ns)), nSamples_stratified(MAX(0,ns_stratified)), samplingType(type) {
	LightToWorld = l2w;
	WorldToLight = l2w.getInverse();
}

const Spectrum Light::LightValue::sample_L(const Vector3<real> &pt, const Vector3<real> &n, real u1, real u2, Vector3<real> *wi, real *pdf, VisibilityTester *vis) const {
	return sample_L(pt, u1, u2, wi, pdf, vis);
}

const real Light::LightValue::pdf(const Vector3<real> &pt, const Vector3<real> &n, const Vector3<real> &wi) const {
	return pdf(pt, wi);
}

const bool VisibilityTester::unOccluded(const Scene *scene) const {
	return !scene->shadowHit(r);
}
