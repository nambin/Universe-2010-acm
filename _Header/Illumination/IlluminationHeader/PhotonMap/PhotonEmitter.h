#ifndef _ILLUMINATION_PHOTON_EMITTER_H_
#define _ILLUMINATION_PHOTON_EMITTER_H_

#include "../../IlluminationSpace.h"
#include "PhotonMap.h"

namespace illuminationspace {

class ILLUMINATIONDLL PhotonEmitter {
public : 
	virtual void emitPhotons(const Scene*) const = 0;
};

class ILLUMINATIONDLL PhotonEmitter_UnbiasedEmitter :public PhotonEmitter {
public : 
	PhotonEmitter_UnbiasedEmitter(size_t numPhoton, PhotonMap* global, PhotonMap* caustic);
	void emitPhotons(const Scene*) const;

	static bool tracePhoton(const Scene*, Ray&, Spectrum& power, PhotonMap *globalPhotonMap, PhotonMap *causticPhotonMap);

public : 
	size_t				numPhotons;

	PhotonMap* const	globalPhotonMap;		// just pointing _ Not Delete
	PhotonMap* const	causticPhotonMap;		// just pointing _ Not Delete
};

};		// illuminationspace

#endif
