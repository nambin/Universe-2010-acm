#pragma once

#include "../../../Principle/Geometry/SRL/SRL3D.h"
#include "IsoSurface_Grid.h"

namespace illuminationspace {

using arithmetic::Grid3;

class ILLUMINATIONDLL IsoSurface_SRLGrid : public IsoSurface_Grid {
protected : 
	class IsoSurface_SRLGridValue;

public:
	IsoSurface_SRLGrid() : IsoSurface_Grid(new IsoSurface_SRLGridValue()) {}
	IsoSurface_SRLGrid(const Transform& o2w, const geometryspace::SRL3D<real>& srl, bool bdrySurface[3][2]) 
	: IsoSurface_Grid(new IsoSurface_SRLGridValue(o2w, srl, bdrySurface)) {}

protected : 
	class ILLUMINATIONDLL IsoSurface_SRLGridValue : public IsoSurface_Grid::IsoSurface_GridValue {
	public : 
		friend class IsoSurface_Grid::IsoSurface_GridValue;

		IsoSurface_SRLGridValue() : IsoSurface_Grid::IsoSurface_GridValue(), interp(true) { init(); }
		IsoSurface_SRLGridValue(const Transform& o2w, const geometryspace::SRL3D<real>&, const bool bdrySurface[3][2]);

		void init();
		void set(const geometryspace::SRL3D<real>&, const bool bdrySurface[3][2]);

		IsoSurface_SRLGridValue * clone() const { return new IsoSurface_SRLGridValue(*this); }

		void getShadingGeometry(const Transform &obj2world, const DifferentialGeometry &dg, DifferentialGeometry &dgShading) const;

		const bool hit(const Ray& r, real& tHit, DifferentialGeometry& geometry) const;
		const bool shadowHit(const Ray& r) const;	

	protected : 
		void setNormals();
		void setSurfaceMarker();

		const bool isSurfaceCell(int i, int j, int k) const { return surfaceMarker(i,j,k);	}		
		const real get_phi(const Vector3<real> &p) const;
		const Vector3<real> get_normal(const Vector3<real>& p) const;

		inline const int PosToVoxel(const Vector3<real> &p, int axis) const;

		inline const real VoxelToPos(int p, int axis) const;
		inline const Vector3<real> VoxelToPos(int x, int y, int z) const;
		inline const Vector3<real> VoxelToPos(const int voxel[3]) const;		

		const bool hitCell(const Ray& ray, const int voxel[3], real t0, real t1, real* tHit, DifferentialGeometry* dg) const;

	private : 
		geometryspace::SRL<real>::SRL_Interpolator				interp;
		geometryspace::SRL<real>::SRL_Lerp						lerp;
		geometryspace::SRL<real>::SRL_Coarse_Lerp				coarse_lerp;

		geometryspace::SRL3D<real>								srl;
		std::vector<geometryspace::SRLNode3D<Vector3<real>>>	normals;

		Grid3<bool,real>										surfaceMarker;
	};
};

inline const int IsoSurface_SRLGrid::IsoSurface_SRLGridValue::PosToVoxel(const Vector3<real> &p, int axis) const {
	const int v = static_cast<int>((p[axis] - bounds.min()[axis]) * InvWidth[axis]);
	return clamp(v, 0, celldim[axis]-1);
}

inline const real IsoSurface_SRLGrid::IsoSurface_SRLGridValue::VoxelToPos(int p, int axis) const {
	return bounds.min()[axis] + p * Width[axis];
}

inline const Vector3<real> IsoSurface_SRLGrid::IsoSurface_SRLGridValue::VoxelToPos(int x, int y, int z) const {
	return bounds.min() + Vector3<real>(x*Width[0], y*Width[1], z*Width[2]);
}

inline const Vector3<real> IsoSurface_SRLGrid::IsoSurface_SRLGridValue::VoxelToPos(const int voxel[3]) const {
	return bounds.min() + Vector3<real>(voxel[0]*Width[0], voxel[1]*Width[1], voxel[2]*Width[2]);
}

};		// illuminationspace