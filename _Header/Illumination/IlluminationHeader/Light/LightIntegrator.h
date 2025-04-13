#ifndef _ILLUMINATION_LIGHT_INTEGRATOR_H_
#define _ILLUMINATION_LIGHT_INTEGRATOR_H_

#include "Light.h"

namespace illuminationspace {

class LightIntegrator {
public : 
	LightIntegrator();
	virtual ~LightIntegrator();

	virtual Spectrum radiance(const Scene *, const Vector3<real>& wo, const Intersection &, const BSDF&) = 0;
	virtual void requestSamples(const Scene*) = 0;
	virtual void getNextSamples();

	static Spectrum directRadiance(const Scene*, const Light&, const Vector3<real>& p, const Vector3<real>& n, const Vector3<real>& wo, const BSDF&
								 , const Sampler&, int lightSample, int numSample);

protected : 
	Sampler*	sampler;
};

class LightIntegrator_UniformSampleAll : public LightIntegrator {
public : 
	LightIntegrator_UniformSampleAll() {}

	Spectrum radiance(const Scene *, const Vector3<real>& wo, const Intersection &, const BSDF&);
	void requestSamples(const Scene*);

protected : 
	Sampler_Offset	offset;
};

class LightIntegrator_UniformSampleOne : public LightIntegrator {
public : 
	LightIntegrator_UniformSampleOne();
	~LightIntegrator_UniformSampleOne();

	Spectrum radiance(const Scene *, const Vector3<real>& wo, const Intersection &, const BSDF&);
	void requestSamples(const Scene*);
};

}

#endif