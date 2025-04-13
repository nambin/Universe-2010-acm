#include <stdafx.h>
#include <Illumination/IlluminationHeader/Integrator/DirectLightingIntegrator.h>
#include <Illumination/IlluminationHeader/Scene/Scene.h>

using namespace illuminationspace;

void DirectLightingIntegrator::requestSamples(const Scene *scene) {	
	SurfaceIntegrator::requestSamples(scene);
}

const Spectrum DirectLightingIntegrator::Li(const Scene *scene, const Ray &ray, real *alpha) const {

	Spectrum retL(0.0f);
	Intersection intersection;

	if(alpha != NULL) 
		*alpha = 1.0f;

	if(!scene->hit(ray, intersection))
		return scene->background;

	Vector3<real> wo = -ray.direction();

	// Area Light
	retL += intersection.Le(wo);

	BSDF *bsdf = intersection.getBSDF(ray);
	if(bsdf == NULL)
		goto returnRadiance;

	// Direct Lighting
	if(!scene->lights.empty() && lightIntegrator) 		
		retL += lightIntegrator->radiance(scene, wo, intersection, *bsdf);

	// Specular Reflection, Transmission
	retL += Li_Specular_Reflection_Transmission(scene, ray, intersection, *bsdf);

returnRadiance :
	Spectrum extinction = ray.extinction();
	if(!extinction.isBlack()) 
		retL *= exp(extinction * intersection.tHit);		

	return retL;
}

