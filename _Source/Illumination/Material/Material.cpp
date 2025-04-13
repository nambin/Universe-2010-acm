#include <stdafx.h>
#include <Illumination/IlluminationHeader/Material/Material.h>

using namespace illuminationspace;

#ifdef USE_OPENMP
MemoryArena BSDF::arena[NUM_OPENMP_THREADS_MAX];
#else
MemoryArena BSDF::arena;
#endif

BSDF *Matte::MatteValue::getBSDF(const DifferentialGeometry& dg) const {
	BSDF *bsdf = BSDF_ALLOC(BSDF)(dg);

	Spectrum r = Kd.evaluate(dg).getClamped();
	real sig = clamp(sigma.evaluate(dg), static_cast<real>(0.0f), static_cast<real>(90.0f));

	if(sig == 0.0f)	bsdf->add(BSDF_ALLOC(Lambertian)(r));
	else			bsdf->add(BSDF_ALLOC(OrenNayar)(r, sig));

	return bsdf;
}

BSDF *Mirror::MirrorValue::getBSDF(const DifferentialGeometry &dg) const {
	BSDF *bsdf = BSDF_ALLOC(BSDF)(dg);

	Spectrum r = Kr.evaluate(dg).getClamped();
	if (!r.isBlack())
		bsdf->add(BSDF_ALLOC(SpecularReflection)(r,	BSDF_ALLOC(FresnelNoOp)()));

	return bsdf;
}

BSDF *Glass::GlassValue::getBSDF(const DifferentialGeometry &dg) const {	
	BSDF *bsdf = BSDF_ALLOC(BSDF)(dg);

	Spectrum r = Kr.evaluate(dg).getClamped();
	Spectrum t = Kt.evaluate(dg).getClamped();

	real ior = index.evaluate(dg);
	if (!r.isBlack()) bsdf->add(BSDF_ALLOC(SpecularReflection)(r, BSDF_ALLOC(FresnelDielectricSchlick)(1.0, ior)));
	if (!t.isBlack()) bsdf->add(BSDF_ALLOC(SpecularTransmission)(t, 1.0, ior, BSDF_ALLOC(FresnelDielectricSchlick)(1.0, ior)));

	//if (!r.isBlack())	bsdf->add(BSDF_ALLOC(SpecularReflection)(r,	BSDF_ALLOC(FresnelDielectric)(1.0, ior)));
	//if (!t.isBlack())	bsdf->add(BSDF_ALLOC(SpecularTransmission)(t, 1.0, ior,	BSDF_ALLOC(FresnelDielectric)(1.0, ior)));

	return bsdf;
}

BSDF *Water::WaterValue::getBSDF(const DifferentialGeometry &dg) const {	
	BSDF *bsdf = BSDF_ALLOC(BSDF)(dg);

	Spectrum r = Kr.evaluate(dg).getClamped();
	Spectrum t = Kt.evaluate(dg).getClamped();

	real ior = index.evaluate(dg);
	if (!r.isBlack()) bsdf->add(BSDF_ALLOC(SpecularReflection)(r, BSDF_ALLOC(FresnelDielectricSchlick)(1.0, ior)));
	if (!t.isBlack()) bsdf->add(BSDF_ALLOC(SpecularTransmission)(t, 1.0, ior, BSDF_ALLOC(FresnelDielectricSchlick)(1.0, ior)));

//	if (!r.isBlack())	bsdf->add(BSDF_ALLOC(SpecularReflection)(r,	BSDF_ALLOC(FresnelDielectric)(1.0, ior)));
//	if (!t.isBlack())	bsdf->add(BSDF_ALLOC(SpecularTransmission)(t, 1.0, ior,	BSDF_ALLOC(FresnelDielectric)(1.0, ior)));

	return bsdf;
}

const Spectrum Water::WaterValue::extinction(const DifferentialGeometry& dg, const Vector3<real>& wi) const {
	Vector3<real> wiLocal(dot(wi, dg.uvw.U), dot(wi, dg.uvw.V), dot(wi, dg.uvw.W));

	if(wiLocal[2] < 0)	return insideExtinction.evaluate(dg);
	else				return outsideExtinction.evaluate(dg);
}

BSDF *Plastic::PlasticValue::getBSDF(const DifferentialGeometry& dg) const {
	BSDF *bsdf = BSDF_ALLOC(BSDF)(dg);

	Spectrum kd = Kd.evaluate(dg).getClamped();
	BxDF *diff = BSDF_ALLOC(Lambertian)(kd);

	Spectrum ks = Ks.evaluate(dg).getClamped();
	Fresnel *fresnel = BSDF_ALLOC(FresnelDielectric)(1.5f, 1.0f);	

	real rough = roughness.evaluate(dg);
	BxDF *spec = BSDF_ALLOC(Microfacet)(ks, fresnel, BSDF_ALLOC(Blinn)(1.0f / rough));

	bsdf->add(diff);
	bsdf->add(spec);
	return bsdf;
}