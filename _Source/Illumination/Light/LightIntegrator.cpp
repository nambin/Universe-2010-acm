#include <stdafx.h>
#include <Illumination/IlluminationHeader/Light/LightIntegrator.h>

using namespace std;
using namespace illuminationspace;

LightIntegrator::LightIntegrator() {
	sampler = new Sampler();
}

LightIntegrator::~LightIntegrator() {
	if(sampler)
		delete sampler;
}

void LightIntegrator::getNextSamples() {
	sampler->getNextSample();
}

Spectrum LightIntegrator::directRadiance(const Scene* scene, const Light& light, const Vector3<real> &p, const Vector3<real> &n, const Vector3<real> &wo, const BSDF &bsdf
									   , const Sampler& samples, int lightSample, int numSample) {
	Spectrum retL(0.0f);

	Vector3<real> wi;
	real pdf;
	VisibilityTester visTester;	

	Spectrum Li = light.sample_L(p, samples.get2D(lightSample,numSample*2), samples.get2D(lightSample,numSample*2+1), &wi, &pdf, &visTester);
	if(pdf > 0.0f && !Li.isBlack()) {
	   Spectrum f = bsdf.f(wo, wi);
	   if(!f.isBlack() && visTester.unOccluded(scene)) {
		   if(light.isDeltaLight()) retL += f * Li * abs(dot(wi, n)) / pdf;
		   else						retL += f * Li * abs(dot(wi, n)) / pdf;
	   }
	}

	return retL;
}
