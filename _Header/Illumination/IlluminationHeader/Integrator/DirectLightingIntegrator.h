#ifndef _ILLUMINATION_DIRECTLIGHTING_INTEGRATOR_H_
#define _ILLUMINATION_DIRECTLIGHTING_INTEGRATOR_H_

#include "SurfaceIntegrator.h"
#include "../Light/LightIntegrator.h"

namespace illuminationspace {

class ILLUMINATIONDLL DirectLightingIntegrator : public SurfaceIntegrator {
public : 
	DirectLightingIntegrator(int _maxDepth=6) 
	: SurfaceIntegrator(_maxDepth) {}

	const Spectrum Li(const Scene*, const Ray&, real* alpha) const;
	void requestSamples(const Scene*);	
};

};		// illuminationspace

#endif