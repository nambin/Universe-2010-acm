#include <stdafx.h>
#include <Illumination/IlluminationHeader/PhotonMap/PhotonGather.h>

using namespace std;
using namespace illuminationspace;

PhotonGather_BasicGather::PhotonGather_BasicGather(size_t num, real _maxDist) {
	numEstimate = num;
	maxDist = _maxDist;
}

const Spectrum PhotonGather_BasicGather::irradiance(const PhotonMap& photonmap, const Vector3<real>& pos, const Vector3<real>& n) const {
	float irrad[3];
	irrad[0] = irrad[1] = irrad[2] = 0;

	NearestPhotons np;
	np.dist2 = (float*)alloca(sizeof(float)*(numEstimate+1));
	np.index = (const Photon**)alloca(sizeof(Photon*)*(numEstimate+1));

	np.pos[0] = pos[0];	np.pos[1] = pos[1];	np.pos[2] = pos[2];
	np.max = static_cast<int>(numEstimate);
	np.found = 0;
	np.got_heap = 0;
	np.dist2[0] = maxDist * maxDist;

	// locate the nearest photons
	photonmap.locate_photons(&np, 1);

	// if less than 8 photons return
	if(np.found < 8)
		return Spectrum::black();

	float pdir[3];
	float normal[3];

	normal[0] = n[0]; normal[1] = n[1]; normal[2] = n[2];

	// sum irradiance from all photons
	for(int i=1;i<=np.found;i++)	{
		const Photon *p = np.index[i];

		// the photon_dir call and following if can be omitted (for speed)
		// if the scene does not have any thin surfaces
		photonmap.photon_dir(pdir, p);
		if((pdir[0]*normal[0] + pdir[1]*normal[1] + pdir[2]*normal[2]) < 0.0f)	{
			irrad[0] += p->power[0];
			irrad[1] += p->power[1];
			irrad[2] += p->power[2];
		}
	}

	const float tmp = (float)((1.0f / PI) / (np.dist2[0]));

	irrad[0] *= tmp;
	irrad[1] *= tmp;
	irrad[2] *= tmp;
	
	return Spectrum(irrad[0], irrad[1], irrad[2]);
}