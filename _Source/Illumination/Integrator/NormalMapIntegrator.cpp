#include <stdafx.h>
#include <Illumination/IlluminationHeader/Integrator/NormalMapIntegrator.h>
#include <Illumination/IlluminationHeader/Primitive/Primitive.h>
#include <Illumination/IlluminationHeader/Scene/Scene.h>

using namespace illuminationspace;

const Spectrum NormalMapIntegrator::Li(const Scene *scene, const Ray &ray, real *alpha) const {
	
	Intersection intersection;

	scene->hit(ray, intersection);

	if(alpha)
		*alpha = 1.0f;

	//return Spectrum(intersection.dg.shadingNormal[0], intersection.dg.shadingNormal[1], intersection.dg.shadingNormal[2]);
	return Spectrum(abs(intersection.dg.shadingNormal[0]), abs(intersection.dg.shadingNormal[1]), abs(intersection.dg.shadingNormal[2]));
}

