#ifndef _ILLUMINATION_MATERIAL_H_
#define _ILLUMINATION_MATERIAL_H_

#include "../Core/Spectrum.h"
#include "../Texture/Texture.h"
#include "../Reflection/BSDF.h"
#include "../Reflection/Reflection.h"
#include "../Reflection/Microfacet.h"

namespace illuminationspace {

class ILLUMINATIONDLL Material {
protected : 
	class MaterialValue;

public :
	Material() : value(NULL) {}	
	Material(MaterialValue *ptr) : value(ptr) {}

	BSDF *getBSDF(const DifferentialGeometry& dg) const { return value->getBSDF(dg); }
	const Spectrum extinction(const DifferentialGeometry& dg, const Vector3<real>& wi) const { return value->extinction(dg, wi); }

protected : 
	class ILLUMINATIONDLL MaterialValue : public basis::RCObject {
	public : 
		MaterialValue() {}
		virtual ~MaterialValue() {}

		virtual MaterialValue *clone() const = 0;	
		virtual BSDF *getBSDF(const DifferentialGeometry& dg) const = 0;
		virtual const Spectrum extinction(const DifferentialGeometry& dg, const Vector3<real>& wi) const { return Spectrum::black(); }
	};

public : 
	//////////////////////////////////////////////////////////////////////////
	// Reference Counted Pointer
	basis::RCPtr<MaterialValue> value;

public : 
	const basis::RCPtr<MaterialValue> &getMaterialValue() const { return value; }	
};

class ILLUMINATIONDLL BlackMaterial : public Material {
protected : 
	class BlackMaterialValue;

public :
	BlackMaterial() : Material(new BlackMaterialValue()) {}

protected : 
	class ILLUMINATIONDLL BlackMaterialValue : public Material::MaterialValue {
	public : 
		BlackMaterialValue * clone() const { return new BlackMaterialValue(*this); }
		BSDF *getBSDF(const DifferentialGeometry& dg) const {
			return NULL;
		}
	};	
}; 

// Matte Class Declarations
class ILLUMINATIONDLL Matte : public Material {
protected : 
	class MatteValue;

public : 
	Matte(Texture<Spectrum> kd, Texture<real> sig) : Material(new MatteValue(kd, sig)) {}

protected : 
	class ILLUMINATIONDLL MatteValue : public Material::MaterialValue {
	public : 
		MatteValue(Texture<Spectrum> kd, Texture<real> sig) { Kd = kd; sigma = sig;	}

		MatteValue * clone() const { return new MatteValue(*this); }
		BSDF *getBSDF(const DifferentialGeometry& dg) const;

	public : 
		Texture<Spectrum>	Kd;
		Texture<real>		sigma;
	};	
};

// Mirror Class Declarations
class ILLUMINATIONDLL Mirror : public Material {
protected : 
	class MirrorValue;

public:
	Mirror(Texture<Spectrum> r) : Material(new MirrorValue(r)) {}

protected : 
	class ILLUMINATIONDLL MirrorValue : public Material::MaterialValue {
	public :
		MirrorValue(Texture<Spectrum> r) { Kr = r;	}

		MirrorValue * clone() const { return new MirrorValue(*this); }
		BSDF *getBSDF(const DifferentialGeometry &dg) const;

	public :
		Texture<Spectrum> Kr;
	};
};

// Glass Class Declarations
class ILLUMINATIONDLL Glass : public Material {
protected : 
	class GlassValue;

public:
	Glass(Texture<Spectrum> r, Texture<Spectrum> t, Texture<real> i) : Material(new GlassValue(r,t,i)) {}

protected : 
	class ILLUMINATIONDLL GlassValue : public Material::MaterialValue {
	public :
		GlassValue(Texture<Spectrum> r, Texture<Spectrum> t, Texture<real> i) {	Kr = r;	Kt = t;	index = i; }

		GlassValue * clone() const { return new GlassValue(*this); }
		BSDF *getBSDF(const DifferentialGeometry &dg) const;

	public :
		Texture<Spectrum>	Kr, Kt;
		Texture<real>		index;
	};
};

// Water Class Declarations
class ILLUMINATIONDLL Water : public Material {
protected : 
	class WaterValue;

public:
	Water(Texture<Spectrum> r, Texture<Spectrum> t, Texture<real> i, Texture<Spectrum> iE, Texture<Spectrum> oE) 
	: Material(new WaterValue(r,t,i,iE,oE)) {}

protected : 
	class ILLUMINATIONDLL WaterValue : public Material::MaterialValue {
	public :
		WaterValue(Texture<Spectrum> r, Texture<Spectrum> t, Texture<real> i, Texture<Spectrum> iE, Texture<Spectrum> oE) 
		{	Kr = r;	Kt = t;	index = i; insideExtinction = iE; outsideExtinction = oE; }

		WaterValue * clone() const { return new WaterValue(*this); }
		BSDF *getBSDF(const DifferentialGeometry &dg) const;
		const Spectrum extinction(const DifferentialGeometry& dg, const Vector3<real>& wi) const;

	public :
		Texture<Spectrum>	Kr, Kt;
		Texture<real>		index;

		Texture<Spectrum>	insideExtinction;
		Texture<Spectrum>	outsideExtinction;
	};
};

// Plastic Class Declarations
class ILLUMINATIONDLL Plastic : public Material {
protected : 
	class PlasticValue;

public :
	Plastic(Texture<Spectrum> kd, Texture<Spectrum> ks, Texture<real> rough)
	: Material(new PlasticValue(kd, ks, rough)) {}

protected : 
	class ILLUMINATIONDLL PlasticValue : public Material::MaterialValue {
	public :
		PlasticValue(Texture<Spectrum> kd, Texture<Spectrum> ks, Texture<real> rough) {	Kd = kd; Ks = ks; roughness = rough; }

		PlasticValue * clone() const { return new PlasticValue(*this); }
		BSDF *getBSDF(const DifferentialGeometry& dg) const;

	public :
		Texture<Spectrum>	Kd, Ks;
		Texture<real>		roughness;
	};
};

};

#endif
