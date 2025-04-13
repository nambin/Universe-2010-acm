#ifndef _ILLUMINATION_FRESNEL_H_
#define _ILLUMINATION_FRESNEL_H_

#include "../Core/Spectrum.h"

namespace illuminationspace {

class ILLUMINATIONDLL Fresnel {
public : 
	virtual const Spectrum evaluate(real cosi) const = 0;

	static const Spectrum FrDiel(real cosi, real cost, const Spectrum &etai, const Spectrum &etat);
	static const Spectrum FrCond(real cosi, const Spectrum &eta, const Spectrum &k);
	static const Spectrum FresnelApproxEta(const Spectrum &Fr);
	static const Spectrum FresnelApproxK(const Spectrum &Fr);
};

class ILLUMINATIONDLL FresnelConductor : public Fresnel {
public :
	FresnelConductor(const Spectrum &e, const Spectrum &kk) : eta(e), k(kk) {};
	const Spectrum evaluate(real cosi) const;

public :
	Spectrum eta, k;
};

class ILLUMINATIONDLL FresnelDielectric : public Fresnel {
public : 
	FresnelDielectric(real ei, real et) : eta_i(ei), eta_t(et) {}
	const Spectrum evaluate(real cosi) const;

public :
	real eta_i, eta_t;
};

// More Think about it _ for symmetry
class ILLUMINATIONDLL FresnelDielectricSchlick : public FresnelDielectric {
public : 
	FresnelDielectricSchlick(real ei, real et) : FresnelDielectric(ei, et) {
		R0 = (eta_t - eta_i) / (eta_t + eta_i);
		R0 *= R0;
	};
	const Spectrum evaluate(real cosi) const;

public :	
	real R0;
};

class ILLUMINATIONDLL FresnelNoOp : public Fresnel {
public : 
	const Spectrum evaluate(real cosi) const { return Spectrum(1.0f); }
};

};		// illuminationspace

#endif