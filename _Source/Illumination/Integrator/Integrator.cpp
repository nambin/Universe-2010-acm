#include <stdafx.h>
#include <Illumination/IlluminationHeader/Integrator/Integrator.h>
#include <Illumination/IlluminationHeader/Scene/Scene.h>
#include <Illumination/IlluminationHeader/Illumination/BasicLibrary.h>

using namespace illuminationspace;

const Spectrum Integrator::Li_Specular_Reflection_Transmission(const Scene* scene, const Ray& ray, const Intersection& isect, const BSDF& bsdf) const {
	Spectrum retL(0.0f);	

	if(ray.depth++ < maxDepth) {
		Vector3<real> wo = - ray.direction(), wi;	
		real pdf; BxDFType type;

		Spectrum f = bsdf.sample_f(wo, &wi, random(), random(), random(), &pdf, BxDFType(BSDF_REFLECTION | BSDF_SPECULAR), &type);
		Spectrum extinction = isect.extinction(wi);

		if(!f.isBlack() && pdf != 0.0f) 
		{
			Ray r(isect.dg.hit_p, wi, extinction, ray.depth);
			r.setTmin(Ray::Epsilon());
			r.setTmax(Ray::Large());

			retL += scene->Li(r) * f * abs(dot(wi, isect.dg.shadingNormal));
		}

		f = bsdf.sample_f(wo, &wi, random(), random(), random(), &pdf, BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR), &type);
		extinction = isect.extinction(wi);

		if(!f.isBlack() && pdf != 0.0f) 
		{
			Ray r(isect.dg.hit_p, wi, extinction, ray.depth);
			r.setTmin(Ray::Epsilon());
			r.setTmax(Ray::Large());

			retL += scene->Li(r) * f * abs(dot(wi, isect.dg.shadingNormal));	
		}
	}

	--ray.depth;

	return retL;
}

//Spectrum Integrator::randomSampleIndirectRadiance(const Scene *scene, const Vector3<real> &wo, const Intersection &isect, const BSDF &bsdf, const Sampler *samples) const {
//	Spectrum retL(0.0f);
//	Vector3<real> wi;
//	real pdf;
//	BxDFType sampledType;
//
//	Spectrum f = bsdf.sample_f(wo, &wi, random(), random(), random(), &pdf, BxDFType(BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_REFLECTION | BSDF_TRANSMISSION), &sampledType);
//
//	if(f.isBlack() || pdf == 0.0f)
//		return Spectrum(0.0f);
//
//	Ray r(isect.dg.hit_p, wi);
//	retL = f * scene->Li(r) * abs(dot(wi, isect.dg.shadingNormal)) / pdf;	
//
//	return retL;
//}
