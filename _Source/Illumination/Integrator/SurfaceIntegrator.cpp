#include <stdafx.h>
#include <Illumination/IlluminationHeader/Integrator/SurfaceIntegrator.h>

using namespace std;
using namespace illuminationspace;

SurfaceIntegrator::SurfaceIntegrator(int _maxDepth) : Integrator(_maxDepth), lightIntegrator(NULL) {}
SurfaceIntegrator::~SurfaceIntegrator() { 
	if(lightIntegrator) 
		delete lightIntegrator; 
}

void SurfaceIntegrator::set_lightIntegrator(LightIntegrator* litg) { 
	if(lightIntegrator)	
		delete lightIntegrator;
	lightIntegrator = litg; 
}

void SurfaceIntegrator::requestSamples(const Scene* scene) {
	if(lightIntegrator)
		lightIntegrator->requestSamples(scene);
}

void SurfaceIntegrator::getNextSamples() { 
	if(lightIntegrator) 
		lightIntegrator->getNextSamples(); 
}

