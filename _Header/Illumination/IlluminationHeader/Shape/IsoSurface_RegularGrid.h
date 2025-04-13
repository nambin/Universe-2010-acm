#ifndef _ILLUMINATION_ISOSURFACE_REGULAR_GRID_H_
#define _ILLUMINATION_ISOSURFACE_REGULAR_GRID_H_

#include "../../../Principle/Arithmetic/GridLibrary/Grid3.h"
#include "IsoSurface_Grid.h"

namespace illuminationspace {

using arithmetic::Grid3;

class ILLUMINATIONDLL IsoSurface_RegularGrid : public IsoSurface_Grid {
protected : 
	class IsoSurface_RegularGridValue;

public:
	IsoSurface_RegularGrid() : IsoSurface_Grid(new IsoSurface_RegularGridValue()) {}
	IsoSurface_RegularGrid(const Transform& o2w, const Grid3<real,real>& grid, bool bdrySurface[3][2]) 
	: IsoSurface_Grid(new IsoSurface_RegularGridValue(o2w, grid, bdrySurface)) {}

protected : 
	class ILLUMINATIONDLL IsoSurface_RegularGridValue : public IsoSurface_Grid::IsoSurface_GridValue {
	public : 
		friend class IsoSurface_Grid::IsoSurface_GridValue;

		IsoSurface_RegularGridValue() : IsoSurface_Grid::IsoSurface_GridValue() { init(); }
		IsoSurface_RegularGridValue(const Transform& o2w, const Grid3<real,real>& phi, const bool bdrySurface[3][2]);

		void init();
		void set(const Grid3<real,real>& phi, const bool bdrySurface[3][2]);

		IsoSurface_RegularGridValue * clone() const { return new IsoSurface_RegularGridValue(*this); }

		void getShadingGeometry(const Transform &obj2world, const DifferentialGeometry &dg, DifferentialGeometry &dgShading) const;

		const bool hit(const Ray& r, real& tHit, DifferentialGeometry& geometry) const;
		const bool shadowHit(const Ray& r) const;	

	protected : 
		void setNormals();
		void setSurfaceMarker();

		const bool isSurfaceCell(int i, int j, int k) const				{ return surfaceMarker(i,j,k);		}

		const real get_phi(const Vector3<real> &p) const				{ return phi(p);					}		
		const Vector3<real> get_normal(const Vector3<real> &p) const	{ return unitVector(normals(p));	}	

		inline const int PosToVoxel(const Vector3<real> &p, int axis) const;

		inline const real VoxelToPos(int p, int axis) const;
		inline const Vector3<real> VoxelToPos(int x, int y, int z) const;
		inline const Vector3<real> VoxelToPos(const int voxel[3]) const;		
		
		const bool hitCell(const Ray& ray, const int voxel[3], real t0, real t1, real* tHit, DifferentialGeometry* dg) const;

	private : 
		Grid3<real,real>			phi;
		Grid3<Vector3<real>,real>	normals;

		Grid3<bool,real>			surfaceMarker;
	};
};

inline const int IsoSurface_RegularGrid::IsoSurface_RegularGridValue::PosToVoxel(const Vector3<real> &p, int axis) const {
	const int v = static_cast<int>((p[axis] - bounds.min()[axis]) * InvWidth[axis]);
	return clamp(v, 0, celldim[axis]-1);
}

inline const real IsoSurface_RegularGrid::IsoSurface_RegularGridValue::VoxelToPos(int p, int axis) const {
	return bounds.min()[axis] + p * Width[axis];
}

inline const Vector3<real> IsoSurface_RegularGrid::IsoSurface_RegularGridValue::VoxelToPos(int x, int y, int z) const {
	return bounds.min() + Vector3<real>(x*Width[0], y*Width[1], z*Width[2]);
}

inline const Vector3<real> IsoSurface_RegularGrid::IsoSurface_RegularGridValue::VoxelToPos(const int voxel[3]) const {
	return bounds.min() + Vector3<real>(voxel[0]*Width[0], voxel[1]*Width[1], voxel[2]*Width[2]);
}

};		// illuminationspace

#endif