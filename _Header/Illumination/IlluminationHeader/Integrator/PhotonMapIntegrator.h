#ifndef _ILLUMINATION_PHOTONMAP_INTEGRATOR_H_
#define _ILLUMINATION_PHOTONMAP_INTEGRATOR_H_

#include "SurfaceIntegrator.h"
#include "../PhotonMap/PhotonEmitter.h"
#include "../PhotonMap/PhotonGather.h"

namespace illuminationspace {

class ILLUMINATIONDLL PhotonMapIntegrator : public SurfaceIntegrator {
public : 
	PhotonMapIntegrator(int _maxDepth=6);

	void init_UnbiasedPhotonEmitter(size_t num);
	void init_BasicPhotonGather(size_t num, real maxDist);
	
	const Spectrum Li(const Scene*, const Ray&, real *alpha) const;
	const Spectrum Li_Photon(const Scene*, const Ray&, real *alpha) const;

	void preProcess(const Scene*);	
	void requestSamples(const Scene*);	

public :
	PhotonMap*		globalPhotonMap;
	PhotonMap*		causticPhotonMap;

	PhotonEmitter*	photonEmitter;
	PhotonGather*	photonGather;

	bool			directPhotonRendering;
};

};		// illuminationspace

#endif