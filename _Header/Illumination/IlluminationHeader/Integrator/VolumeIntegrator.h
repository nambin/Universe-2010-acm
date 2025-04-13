#ifndef _ILLUMINATION_VOLUME_INTEGRATOR_H_
#define _ILLUMINATION_VOLUME_INTEGRATOR_H_

#include "Integrator.h"
#include "../Light/LightIntegrator.h"

namespace illuminationspace {

class ILLUMINATIONDLL VolumeIntegrator : public Integrator {
public : 
	VolumeIntegrator(int _maxDepth) : Integrator(_maxDepth) {}
};

};		// illuminationspace

#endif