#ifndef _ILLUMINATION_PHOTON_GATHER_H_
#define _ILLUMINATION_PHOTON_GATHER_H_

#include "../../IlluminationSpace.h"
#include "PhotonMap.h"
#include "../Core/Spectrum.h"

namespace illuminationspace {

class ILLUMINATIONDLL PhotonGather {
public : 
	virtual const Spectrum irradiance(const PhotonMap&, const Vector3<real>& pos, const Vector3<real>& n) const = 0;
};

class ILLUMINATIONDLL PhotonGather_BasicGather : public PhotonGather {
public : 
	PhotonGather_BasicGather(size_t num, real maxDist);
	const Spectrum irradiance(const PhotonMap&, const Vector3<real>& pos, const Vector3<real>& n) const;

public : 
	size_t		numEstimate;
	real		maxDist;
};

};			// illuminationspace

#endif