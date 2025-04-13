#ifndef _ILLUMINATION_BSDF_H_
#define _ILLUMINATION_BSDF_H_

#include "../Illumination/MemoryArena.h"
#include "../Shape/Shape.h"
#include "Reflection.h"

namespace illuminationspace {

#define MAX_BxDFS 8
#define BSDF_ALLOC(T) new (BSDF::Alloc(sizeof(T))) T

class BxDF;
struct DifferentialGeometry;

class BSDF {
public : 
	inline BSDF(const DifferentialGeometry& _dg, real _eta = 1.0f);

private : 
	~BSDF() {}

public : 
	inline void add(BxDF *b);
	inline const int numComponents() const { return nBxDFs; }
	inline const int numComponents(BxDFType flags) const;

	inline const Vector3<real> WorldToLocal(const Vector3<real> &_vv) const;
	inline const Vector3<real> LocalToWorld(const Vector3<real> &_vv) const;

	const Spectrum f(const Vector3<real> &woW, const Vector3<real> &wiW, BxDFType flags = BSDF_ALL) const; 
	const Spectrum sample_f(const Vector3<real> &woW, Vector3<real> *wiW, real u1, real u2, real u3, real *pdf, BxDFType flags, BxDFType *sampledType) const;
	// Check Routine's validity
	
	const bool reflectPhoton(const Vector3<real> &wIn, Vector3<real> *wOut, real u1, real u2, real u3, const Spectrum &power, Spectrum *retPower, BxDFType *reflectType);
	
	const Spectrum rho(BxDFType flags = BSDF_ALL) const;
	const Spectrum rho(const Vector3<real> &wo, BxDFType flags = BSDF_ALL) const;

public : 
	const DifferentialGeometry	dg;
	const real					eta;

	Vector3<real>				u, v, n, nShading;

	int							nBxDFs;
	BxDF*						bxdfs[MAX_BxDFS];

#ifdef USE_OPENMP
	static void *Alloc(unsigned int sz) { return arena[omp_get_thread_num()].Alloc(sz); }
	static void FreeAll()				{ return arena[omp_get_thread_num()].FreeAll(); }

	static MemoryArena					arena[NUM_OPENMP_THREADS_MAX];	
#else
	static void *Alloc(unsigned int sz) { return arena.Alloc(sz); }
	static void FreeAll()				{ return arena.FreeAll(); }

	static MemoryArena					arena;	
#endif
};

inline BSDF::BSDF(const DifferentialGeometry& _dg, real _eta) 
: dg(_dg), eta(_eta) {
	nBxDFs = 0;
	u = dg.uvw.U;
	v = dg.uvw.V;
	n = dg.uvw.W;
	nShading = dg.shadingNormal;
}

inline void BSDF::add(BxDF *b) {
	bxdfs[nBxDFs++] = b;
}

inline const int BSDF::numComponents(BxDFType flags) const {
	int num = 0;
	for (int i = 0; i < nBxDFs; ++i)
		if (bxdfs[i]->MatchesFlags(flags)) ++num;
	return num;
}

inline const Vector3<real> BSDF::WorldToLocal(const Vector3<real> &_vv) const {
	return Vector3<real>(dot(_vv, u), dot(_vv, v), dot(_vv, n));
}
inline const Vector3<real> BSDF::LocalToWorld(const Vector3<real> &_vv) const {
	return Vector3<real>(u.x() * _vv.x() + v.x() * _vv.y() + n.x() * _vv.z(), 
						 u.y() * _vv.x() + v.y() * _vv.y() + n.y() * _vv.z(),
						 u.z() * _vv.x() + v.z() * _vv.y() + n.z() * _vv.z());
}

};		// illuminationspace

#endif