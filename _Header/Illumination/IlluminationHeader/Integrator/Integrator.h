#ifndef _ILLUMINATION_INTEGRATOR_H_
#define _ILLUMINATION_INTEGRATOR_H_

#include "../Core/Sampler.h"
#include "../Core/Spectrum.h"
#include "../Core/Ray.h"
#include "../Light/Light.h"
#include "../Material/Material.h"
#include "../Primitive/Primitive.h"
#include "../PhotonMap/PhotonMap.h"

namespace illuminationspace {

class Scene;
class ILLUMINATIONDLL Integrator {
public :
	Integrator(int _maxDepth) : maxDepth(_maxDepth) {}

	virtual const Spectrum Li(const Scene *scene, const Ray &ray, real *alpha) const = 0;
	const Spectrum Li_Specular_Reflection_Transmission(const Scene*, const Ray&, const Intersection&, const BSDF&) const;

	virtual void requestSamples(const Scene *scene) = 0;
	virtual void getNextSamples() = 0;

	virtual void preProcess(const Scene *scene) { requestSamples(scene); }
	virtual void printStatus(const Scene *scene) {}
	virtual void postProcess(const Scene *scene) {}
		
//	Spectrum randomSampleIndirectRadiance(const Scene *, const Vector3<real> &, const Intersection &, const BSDF &, const Sampler *sample) const;

protected : 
	const int	maxDepth;
};

};		// illuminationspace

#endif