#include <stdafx.h>
#include <Illumination/IlluminationHeader/Light/LightIntegrator.h>

using namespace std;
using namespace illuminationspace;

Spectrum LightIntegrator_UniformSampleAll::radiance(const Scene* scene, const Vector3<real>& wo, const Intersection& isect, const BSDF& bsdf) {
	Spectrum retL(0.0f);

	for(unsigned int i=0;i<scene->lights.size();++i) {
		const Light& light = scene->lights[i];
		const unsigned int nSamples = sampler->num2D(offset[i]);

		Spectrum L_light(0.0f);

		for(unsigned int j=0;j<nSamples;++j)
			L_light += directRadiance(scene, light, isect.dg.hit_p, isect.dg.shadingNormal, wo, bsdf, *sampler, offset[i], j);

		retL += L_light / (real)nSamples;
	}

	return retL;
}

void LightIntegrator_UniformSampleAll::requestSamples(const Scene* scene) {
	offset.resize(scene->lights.size());
	for(size_t i=0;i<scene->lights.size();++i) 
		offset.set( i , sampler->add2D(scene->lights[i].nSamples(), scene->lights[i].nSamples_stratified(), scene->lights[i].samplingType()) );

	sampler->allocate();
}