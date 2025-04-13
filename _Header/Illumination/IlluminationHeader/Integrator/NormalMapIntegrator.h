#ifndef _NORMAL_MAP_INTEGRATOR_H_
#define _NORMAL_MAP_INTEGRATOR_H_

#include "SurfaceIntegrator.h"

namespace illuminationspace {

class ILLUMINATIONDLL NormalMapIntegrator : public SurfaceIntegrator {
public : 
	NormalMapIntegrator(int _maxDepth = 3) : SurfaceIntegrator(_maxDepth) {}

	const Spectrum Li(const Scene *scene, const Ray &ray, real *alpha) const;
	void requestSamples(const Scene *scene) {}
	void getNextSamples() {}

public : 
};

};		// illuminationspace

#endif
