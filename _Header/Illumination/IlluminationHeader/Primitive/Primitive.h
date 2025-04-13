#ifndef _ILLUMINATION_PRIMITIVE_H_
#define _ILLUMINATION_PRIMITIVE_H_

#include "../Core/Ray.h"
#include "../Shape/Shape.h"
#include "../Material/Material.h"
#include "../Light/AreaLight.h"

namespace illuminationspace {

struct Intersection;

class ILLUMINATIONDLL Primitive {
protected : 
	class PrimitiveValue;

public : 
	Primitive() : value(NULL) {}
	Primitive(PrimitiveValue *ptr) : value(ptr) {}

	const BBox worldBound() const							{ return value->worldBound(); }

	const bool canIntersect() const							{ return value->canIntersect(); }
	void refine(std::vector<Primitive> &refined) const		{ value->refine(refined);		}
	void fullyRefine(std::vector<Primitive> &refined) const	{ value->fullyRefine(refined);	}

	const AreaLight& getAreaLight() const					{ return value->getAreaLight(); }
	BSDF *getBSDF(const DifferentialGeometry &dg) const		{ return value->getBSDF(dg); }

	const bool hit(const Ray& r, Intersection &intersection) const	{ return value->hit(r, intersection);	}
	const bool shadowHit(const Ray& r) const						{ return value->shadowHit(r);			}

protected :
	class ILLUMINATIONDLL PrimitiveValue : public basis::RCObject {
	public : 
		PrimitiveValue() {}

		virtual PrimitiveValue* clone() const;

		virtual const BBox worldBound() const = 0;

		virtual const bool canIntersect() const = 0;
		virtual void refine(std::vector<Primitive> &refined) const;			// called by fullyRefine() function
		void fullyRefine(std::vector<Primitive> &refined) const;			// Used for aggregates

		virtual const AreaLight& getAreaLight() const;
		virtual BSDF* getBSDF(const DifferentialGeometry &dg) const;
		virtual const Spectrum extinction(const DifferentialGeometry& dg, const Vector3<real>& wi) const;

		virtual const bool hit(const Ray& r, Intersection &intersection) const;
		virtual const bool shadowHit(const Ray& r) const;
	};

	//////////////////////////////////////////////////////////////////////////
	// Reference Counted Pointer
	basis::RCPtr<PrimitiveValue> value;

public : 
	const basis::RCPtr<PrimitiveValue>& getPrimitiveValue() const { return value; }	
	friend struct Intersection;
};

struct Intersection	{
public : 
	real tHit;
	DifferentialGeometry dg;
	const Primitive::PrimitiveValue *primitiveValue;

public : 
	const Spectrum Le(const Vector3<real> &wi) const {
		const AreaLight& area = primitiveValue->getAreaLight();
//		if(area.getLightValue())
//			std::cout << "";
		return area.getLightValue() ? area.L(dg.hit_p, dg.uvw.w(), wi) : Spectrum(0.0f);
	}
	BSDF *getBSDF(const Ray &ray) const {		
		return primitiveValue->getBSDF(dg);
	}
	const Spectrum extinction(const Vector3<real>& wi) const {
		return primitiveValue->extinction(dg, wi);
	}
};

};		// illuminationspace

#endif