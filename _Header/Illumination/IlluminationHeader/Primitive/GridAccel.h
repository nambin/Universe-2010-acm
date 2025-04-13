#ifndef _ILLUMINATION_GRID_ACCEL_H_
#define _ILLUMINATION_GRID_ACCEL_H_

#include "../Illumination/ObjectArena.h"
#include "Aggregate.h"
#include "../Illumination/BasicLibrary.h"

namespace illuminationspace {

struct Voxel;
class ILLUMINATIONDLL GridAccel : public Aggregate {
protected :
	class GridAccelValue;

public : 
	GridAccel(const std::vector<Primitive> &p, bool forRefined = false, bool refineImmediately = true)
	: Aggregate(new GridAccelValue(p, forRefined, refineImmediately)) {}

protected : 
	class ILLUMINATIONDLL GridAccelValue : public Aggregate::AggregateValue {
	public : 
		GridAccelValue(const std::vector<Primitive> &p, bool forRefined = false, bool refineImmediately = true);

		const BBox worldBound() const { return bounds; }
		const bool canIntersect() const { return true; }

		const bool hit(const Ray& r, Intersection &intersection) const;
		const bool shadowHit(const Ray& r) const;

		const int PosToVoxel(const Vector3<real> &p, int axis) const {
			//int v = Float2Int((p[axis] - bounds.min()[axis]) * InvWidth[axis]);
			int v = static_cast<int>((p[axis] - bounds.min()[axis]) * InvWidth[axis]);
			return clamp(v, 0, NVoxels[axis]-1);
		}

		const real VoxelToPos(int p, int axis) const {
			return bounds.min()[axis] + p * Width[axis];
		}

		const Vector3<real> VoxelToPos(int x, int y, int z) const {
			return bounds.min() + Vector3<real>(x*Width[0], y*Width[1], z*Width[2]);
		}

		const int Offset(int x, int y, int z) const {
			return x * NVoxels[1] * NVoxels[2] + y * NVoxels[2] + z;
		}

	public : 
		bool gridForRefined;
		unsigned int nMailboxes;
		MailboxPrim *mailboxes;

		BBox bounds;
		int NVoxels[3];

		Vector3<real> Width, InvWidth;
		Voxel **voxels;

		ObjectArena<Voxel> voxelArena;

		mutable int curMailboxId;
	};
};	

struct Voxel {
public : 
	Voxel(MailboxPrim *op) : allCanIntersect(false), nPrimitives(1), onePrimitive(op) {}
	void AddPrimitive(MailboxPrim *prim);

	const bool hit(const Ray& r, Intersection &intersection, int rayId);
	const bool shadowHit(const Ray& r, int rayId);

public : 
	union {
		MailboxPrim *onePrimitive;
		MailboxPrim **primitives;
	};
	unsigned int allCanIntersect:1;
	unsigned int nPrimitives:31;
};

};		// illuminationspace

#endif
