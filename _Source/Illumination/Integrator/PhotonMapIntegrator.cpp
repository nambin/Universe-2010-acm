#include <stdafx.h>
#include <Illumination/IlluminationHeader/Integrator/PhotonMapIntegrator.h>

using namespace illuminationspace;

PhotonMapIntegrator::PhotonMapIntegrator(int _maxDepth) 
: SurfaceIntegrator(_maxDepth), globalPhotonMap(NULL), causticPhotonMap(NULL), photonEmitter(NULL), photonGather(NULL), directPhotonRendering(false) {}

void PhotonMapIntegrator::preProcess(const Scene *scene) {
	SurfaceIntegrator::preProcess(scene);
	photonEmitter->emitPhotons(scene);
}

void PhotonMapIntegrator::requestSamples(const Scene *scene) {
	SurfaceIntegrator::requestSamples(scene);
}

const Spectrum PhotonMapIntegrator::Li(const Scene* scene, const Ray& ray, real* alpha) const {

	if(directPhotonRendering)
		return Li_Photon(scene, ray, alpha);

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

	// Render Caustic Using Caustic PhotonMap
	if(bsdf->numComponents(BxDFType(BSDF_DIFFUSE | BSDF_REFLECTION | BSDF_TRANSMISSION))) {
		Spectrum brdf = bsdf->f(wo, wo, BxDFType(BSDF_DIFFUSE | BSDF_REFLECTION | BSDF_TRANSMISSION));
		Spectrum irrad = photonGather->irradiance(*causticPhotonMap, intersection.dg.hit_p, intersection.dg.uvw.W);

		retL += (irrad * brdf);
	}

	// Specular Reflection, Transmission
	retL += Li_Specular_Reflection_Transmission(scene, ray, intersection, *bsdf);

returnRadiance :
	Spectrum extinction = ray.extinction();
	if(!extinction.isBlack()) 
		retL *= exp(extinction * intersection.tHit);		

	return retL;
}
