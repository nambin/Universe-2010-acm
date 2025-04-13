#ifndef _ILLUMINATION_SURFACE_INTEGRATOR_H_
#define _ILLUMINATION_SURFACE_INTEGRATOR_H_

#include "Integrator.h"
#include "../Light/LightIntegrator.h"

namespace illuminationspace {

class ILLUMINATIONDLL SurfaceIntegrator	: public Integrator {
public :
	SurfaceIntegrator(int _maxDepth);
	virtual ~SurfaceIntegrator();

	void set_lightIntegrator(LightIntegrator*);

	void requestSamples(const Scene*);
	void getNextSamples();

	void preProcess(const Scene* scene) { Integrator::preProcess(scene); }

protected : 
	LightIntegrator*	lightIntegrator;
};

};		// illuminationspace

#endif