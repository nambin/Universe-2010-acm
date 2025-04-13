#ifndef _ILLUMINATION_PHOTON_MAP_H_
#define _ILLUMINATION_PHOTON_MAP_H_

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "../Illumination/Environments.h"
#include "../Shape/Shape.h"

namespace illuminationspace {

#define swapPhotonPtr(ph,a,b)	{ Photon *ph2 = ph[a]; ph[a] = ph[b]; ph[b] = ph2;	}

// This is the photon
// The power is not compressed so the
// size is 28 bytes
typedef struct Photon {
	float			pos[3];         // photon position
	short			plane;          // splitting plane for kd-tree
	unsigned char	theta, phi;		// incoming direction
	float			power[3];       // photon power (uncompressed)
} Photon;

// This structure is used only to locate the
// nearest photons
typedef struct NearestPhotons {
	int				max;
	int				found;
	int				got_heap;
	float			pos[3];
	float *			dist2;
	const Photon**	index;
} NearestPhotons;

class PhotonMap {
public : 
	PhotonMap(int max_phot);
	~PhotonMap();

	void store(const float power[3], const float pos[3], const float dir[3]);
	void scale_photon_power(const float scale);
	void balance(void);

	void locate_photons(NearestPhotons* np, const real pos[3], const real max_dist, const int nPhotons) const;
	void locate_photons(NearestPhotons* const np, const int index) const;
	void photon_dir(float *dir, const Photon *p) const;

	int nStoredPhotons() const { return stored_photons;	}
	const Photon* getPhotons() const { return photons;	}

private :
	void balance_segment(Photon **pbal, Photon **porg, const int index, const int start, const int end);
	void median_split(Photon **p, const int start, const int end, const int median, const int axis);

	Photon*		photons;
	int			stored_photons;
	int			half_stored_photons;
	int			max_photons;
	int			prev_scale;

	float		costheta[256];
	float		sintheta[256];
	float		cosphi[256];
	float		sinphi[256];

	float		bbox_min[3];
	float		bbox_max[3];
};

};			// illuminationspace

#endif