#pragma once

#include "../../../Principle/Arithmetic/GridLibrary/Grid3.h"
#include "IsoSurface.h"

namespace illuminationspace {

using arithmetic::Grid3;

class ILLUMINATIONDLL IsoSurface_Grid : public IsoSurface {
protected : 
	class IsoSurface_GridValue;

public:
	IsoSurface_Grid(IsoSurface_GridValue* ptr) : IsoSurface(ptr) {}	

protected : 
	class ILLUMINATIONDLL IsoSurface_GridValue : public IsoSurface::IsoSurfaceValue {
	public : 
		IsoSurface_GridValue() : IsoSurface::IsoSurfaceValue() { init(); }
		IsoSurface_GridValue(const Transform& o2w) : IsoSurface::IsoSurfaceValue(o2w) { init(); }

		void init();
		void set(const int dim[3], const BBox bounds, const bool bdrySurface[3][2]);

		const BBox objectBound() const { return bounds; }

	protected : 
		template<typename ISO> const bool hitCore(const ISO&, const Ray&, real* tHit=0, DifferentialGeometry* dg=0) const;
		template<typename ISO> const bool hitBoundarySurface(const ISO&, const Ray&, real rayT, const int voxel[3], real* tHit, DifferentialGeometry* dg) const;
		template<typename ISO> const bool hitBoundarySurface(const ISO&, const Ray&, real rayT, const int voxel[3], int axis, real* tHit, DifferentialGeometry* dg) const;

		virtual const bool hitCell(const Ray &ray, const int voxel[3], real t0, real t1, real *tHit, DifferentialGeometry *dg) const = 0;

	protected : 
		int				dim[3];				// number of grid points (samples) in each direction
		int				celldim[3];			// number of grid cells (one less than samples)
		BBox			bounds;				// object-space extent 

		Vector3<real>	Width, InvWidth;

		bool			isBoundarySurface[3][2];
	};
};

};		// illuminationspace

#include "IsoSurface_Grid_Intersect.h"