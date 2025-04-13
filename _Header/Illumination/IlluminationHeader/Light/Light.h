#ifndef _ILLUMINATION_LIGHT_H_
#define _ILLUMINATION_LIGHT_H_

#include "../../../Basis/MemoryManagement/RCPtr.h"
#include "../Illumination/Random.h"
#include "../Core/Transform.h"
#include "../Core/Ray.h"
#include "../Core/Sampler.h"

namespace illuminationspace {

class Scene;
struct VisibilityTester;

class ILLUMINATIONDLL Light {
protected :
	class LightValue;

public :
	Light() : value(NULL) {}
	Light(LightValue *ptr) : value(ptr) {}

	const Spectrum power(const Scene *scene) const	{ return value->power(scene);	}
	const bool isDeltaLight() const					{ return value->isDeltaLight();	}

	const real pdf(const Vector3<real> &pt, const Vector3<real> &wi) const								// Need to check _ When is used ?
	{ return value->pdf(pt, wi);	}

	const real pdf(const Vector3<real> &pt, const Vector3<real> &n, const Vector3<real> &wi) const		// Not used
	{ return value->pdf(pt, n, wi);	}	

	const Spectrum sample_L(const Vector3<real> &pt, Vector3<real> *wi, VisibilityTester *vis) const 
	{ return value->sample_L(pt, wi, vis);	}

	const Spectrum sample_L(const Vector3<real> &pt, real u1, real u2, Vector3<real> *wi, real *pdf, VisibilityTester *vis) const 
	{ return value->sample_L(pt, u1, u2, wi, pdf, vis);	}

	const Spectrum sample_L(const Vector3<real> &pt, const Vector3<real> &n, real u1, real u2, Vector3<real> *wi, real *pdf, VisibilityTester *vis) const
	{ return value->sample_L(pt, n, u1, u2, wi, pdf, vis);	}											// Not used

	void generatePhotonPath(Vector3<real> &pos, Vector3<real> &dir) const 
	{ return value->generatePhotonPath(pos, dir);	}

	//////////////////////////////////////////////////////////////////////////
	// Getter
	const Transform LightToWorld() const	{ return value->LightToWorld; }
	const Transform WorldToLight() const	{ return value->WorldToLight; }
	const int nSamples() const				{ return value->nSamples; }
	const int nSamples_stratified() const	{ return value->nSamples_stratified; }
	const Sampler::SAMPLING_TYPE samplingType() const { return value->samplingType; }

protected : 
	class ILLUMINATIONDLL LightValue : public basis::RCObject {
	public : 
		LightValue();
		LightValue(const Transform &l2w, int ns=1, int ns_stratified=0, Sampler::SAMPLING_TYPE type=Sampler::RANDOM);
		virtual ~LightValue() {}

		virtual LightValue *clone() const = 0;		

		virtual const Spectrum power(const Scene *) const = 0;
		virtual const bool isDeltaLight() const = 0;

		virtual const real pdf(const Vector3<real> &pt, const Vector3<real> &wi) const = 0;								// Need to check _ When is used ?
		virtual const real pdf(const Vector3<real> &pt, const Vector3<real> &n, const Vector3<real> &wi) const;			// Not used

		virtual const Spectrum sample_L(const Vector3<real> &pt, Vector3<real> *wi, VisibilityTester *vis) const = 0;
		virtual const Spectrum sample_L(const Vector3<real> &pt, real u1, real u2, Vector3<real> *wi, real *pdf, VisibilityTester *vis) const = 0;
		virtual const Spectrum sample_L(const Vector3<real> &pt, const Vector3<real> &n, real u1, real u2, Vector3<real> *wi, real *pdf, VisibilityTester *vis) const;	
		// Not used

		virtual void generatePhotonPath(Vector3<real> &pos, Vector3<real> &dir) const = 0;

	public :
		Transform				LightToWorld, WorldToLight;

		const int				nSamples;
		const int				nSamples_stratified;
		Sampler::SAMPLING_TYPE	samplingType;	
	};

	//////////////////////////////////////////////////////////////////////////
	// Reference Counted Pointer
	basis::RCPtr<LightValue> value;

public : 
	const basis::RCPtr<LightValue> &getLightValue() const { return value; }
};

struct ILLUMINATIONDLL VisibilityTester {
public : 
	void setSegment(const Vector3<real> &p1, const Vector3<real> &p2)	{
		r.setOrigin(p1);
		r.setDirection(p2-p1);
		r.setTmin(Ray::Epsilon());
		r.setTmax(1.0f - Ray::Epsilon());		
	}
	void setRay(const Vector3<real> &p, const Vector3<real> &w)	{
		r.setOrigin(p);
		r.setDirection(w);
		r.setTmin(Ray::Epsilon());
		r.setTmax(Ray::Large());
	}
	const bool unOccluded(const Scene *scene) const;

public : 
	Ray r;
};

};		// illuminationspace

#endif