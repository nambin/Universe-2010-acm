#include <stdafx.h>
#include <Illumination/IlluminationHeader/Integrator/PhotonMapIntegrator.h>

using namespace std;
using namespace illuminationspace;

void PhotonMapIntegrator::init_UnbiasedPhotonEmitter(size_t numPhotons) {
	int maxStoredPhotons = numPhotons * 2.5f;

	globalPhotonMap = new PhotonMap(maxStoredPhotons);
	causticPhotonMap = new PhotonMap(maxStoredPhotons);

	photonEmitter = new PhotonEmitter_UnbiasedEmitter(numPhotons, globalPhotonMap, causticPhotonMap);	
}

void PhotonMapIntegrator::init_BasicPhotonGather(size_t num, real maxDist) {
	photonGather = new PhotonGather_BasicGather(num, maxDist);
}