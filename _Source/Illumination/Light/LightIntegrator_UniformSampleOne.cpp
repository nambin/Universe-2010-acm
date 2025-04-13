#include <stdafx.h>
#include <Illumination/IlluminationHeader/Light/LightIntegrator.h>

using namespace std;
using namespace illuminationspace;

LightIntegrator_UniformSampleOne::LightIntegrator_UniformSampleOne() {}
LightIntegrator_UniformSampleOne::~LightIntegrator_UniformSampleOne() {}

Spectrum LightIntegrator_UniformSampleOne::radiance(const Scene* scene, const Vector3<real>& wo, const Intersection& isect, const BSDF& bsdf) {
	int nLights = (int)scene->lights.size(), lightNum;

	lightNum = Floor2Int(sampler->get1D(0,0) * nLights);	
	lightNum = MIN(lightNum, nLights-1);

	const Light& light = scene->lights[lightNum];
	const unsigned int nSamples = sampler->num2D(0);

	Spectrum retL(0.0f), L_light(0.0f);

	L_light += directRadiance(scene, light, isect.dg.hit_p, isect.dg.shadingNormal, wo, bsdf, *sampler, 0, 0);
	retL += L_light / (real)nSamples;

	return nLights * retL;
}

void LightIntegrator_UniformSampleOne::requestSamples(const Scene* scene) {	
	sampler->add1D(1);
	sampler->add2D(1);
	sampler->allocate();
}