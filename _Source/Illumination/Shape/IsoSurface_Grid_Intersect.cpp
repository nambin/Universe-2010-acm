#include <stdafx.h>
#include <Illumination/IlluminationHeader/Shape/IsoSurface_Grid.h>

using namespace illuminationspace;
//
//const bool IsoSurface_Grid::IsoSurface_GridValue::hitCore(const Ray &r, real *tHit, DifferentialGeometry *dg) const {
//	real rayT, tmax;
//
//	if (bounds.isInside(r(r.tMin())))			rayT = r.tMin();
//	else if (!bounds.intersectP(r, rayT, tmax))	return false;
//
//	Vector3<real> gridIntersect = r(rayT);
//
//	///////////////////////////////////////////////////////////////////////////////
//	// Set up 3D DDA for ray
//	real NextCrossingT[3], DeltaT[3];
//	int Step[3], Out[3], Voxel[3];
//
//	for(int axis=0; axis<3; ++axis){
//
//		///////////////////////////////////////////////////////////////////////////////
//		// Compute current voxel for axis
//		Voxel[axis] = PosToVoxel(gridIntersect, axis);
//
//		if (r.direction()[axis] >= 0) {
//			///////////////////////////////////////////////////////////////////////////////
//			// Handle ray with Positive direction for voxel stepping
//			NextCrossingT[axis] = rayT + (VoxelToPos(Voxel[axis]+1, axis) - gridIntersect[axis]) * r.invDirection()[axis];
//			DeltaT[axis] = Width[axis] * r.invDirection()[axis];
//			Step[axis] = 1;
//			Out[axis] = celldim[axis];
//		}
//		else {
//			///////////////////////////////////////////////////////////////////////////////
//			// Handle ray with negative direction for voxel stepping
//			NextCrossingT[axis] = rayT + (VoxelToPos(Voxel[axis], axis) - gridIntersect[axis]) * r.invDirection()[axis];
//			DeltaT[axis] = -Width[axis] * r.invDirection()[axis];
//			Step[axis] = -1;
//			Out[axis] = -1;
//		}
//	}
//
//	if(hitBoundarySurface(r, rayT, Voxel, tHit, dg))
//		return true;
//
//	///////////////////////////////////////////////////////////////////////////////
//	// Walk ray through grid
//	while(true) {
//		///////////////////////////////////////////////////////////////////////////////
//		// Find _stepAxis_ for stepping to next voxel
//		int bits =	((NextCrossingT[0] < NextCrossingT[1]) << 2) +
//			((NextCrossingT[0] < NextCrossingT[2]) << 1) +
//			((NextCrossingT[1] < NextCrossingT[2]));
//
//		const int cmpToAxis[8] = { 2, 1, 2, 1, 2, 2, 0, 0 };
//
//		int stepAxis = cmpToAxis[bits];
//
//		// check this cell
//		if(surfaceMarker(Voxel[0], Voxel[1], Voxel[2]) && hitCell(r, Voxel, rayT, MIN(r.tMax(), NextCrossingT[stepAxis]), tHit, dg))
//			return true;
//
//		if(r.tMax() < NextCrossingT[stepAxis])
//			break;
//
//		Voxel[stepAxis] += Step[stepAxis];
//		if(Voxel[stepAxis] == Out[stepAxis]) {
//			Voxel[stepAxis] -= Step[stepAxis];
//			if(hitBoundarySurface(r, NextCrossingT[stepAxis], Voxel, stepAxis, tHit, dg))
//				return true;
//			break;
//		}
//
//		rayT = NextCrossingT[stepAxis];
//		NextCrossingT[stepAxis] += DeltaT[stepAxis];
//	}
//
//	return false;
//}
//
//const bool IsoSurface_Grid::IsoSurface_GridValue::hitBoundarySurface(const Ray &ray, real rayT, const int voxel[3], real *tHit, DifferentialGeometry *dg) const {
//	for(int axis=0;axis<3;++axis) {
//		if(hitBoundarySurface(ray, rayT, voxel, axis, tHit, dg)) {
//			return true;
//		}
//	}
//
//	return false;
//}
//
//const bool IsoSurface_Grid::IsoSurface_GridValue::hitBoundarySurface(const Ray &ray, real rayT, const int voxel[3], int axis, real *tHit, DifferentialGeometry *dg) const {
//
//	static const real eps = 1.e-4;
//
//	int flagBdryVoxel;	
//	if(voxel[axis] == 0 && isBoundarySurface[axis][0])
//		flagBdryVoxel = 0;
//	else if(voxel[axis] == celldim[axis]-1 && isBoundarySurface[axis][1]) 
//		flagBdryVoxel = 1;
//	else
//		return false;
//
//	Vector3<real> pos = ray(rayT);
//
//	if(abs(pos[axis] - bounds.pp[flagBdryVoxel][axis]) > eps * Width[axis])
//		return false;
//
//	bounds.constrainedBdry(pos);
//	if(val(pos) < 0) {
//		if(!tHit)
//			return true;
//
//		*tHit = rayT;
//
//		dg->hit_p = ObjectToWorld.transformPoint(pos);
//		dg->shadingNormal = Vector3<real>(0.0f);
//		dg->shadingNormal[axis] = flagBdryVoxel * 2.0f - 1.0f;
//		dg->shadingNormal = ObjectToWorld.transformNormal(dg->shadingNormal);
//		dg->uvw.initFromUnitVector(dg->shadingNormal);
//
//		return true;
//	}
//
//	return false;
//}
