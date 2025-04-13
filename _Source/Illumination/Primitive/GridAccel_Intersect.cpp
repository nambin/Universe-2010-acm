#include <stdafx.h>
#include <Illumination/IlluminationHeader/Primitive/GridAccel.h>

using namespace illuminationspace;

const bool GridAccel::GridAccelValue::hit(const Ray& r, Intersection &intersection) const {

	///////////////////////////////////////////////////////////////////////////////
	// Check ray against overall grid bounds
	real rayT, tmax;

	if (bounds.isInside(r(r.tMin())))			rayT = r.tMin();
	else if (!bounds.intersectP(r, rayT, tmax))	return false;

	Vector3<real> gridIntersect = r(rayT);

	///////////////////////////////////////////////////////////////////////////////
	// Get ray mailbox id
//	int rayId = ++curMailboxId;
#ifdef USE_OPENMP
	int rayId;
{
#pragma omp critical 
	{
		rayId = curMailboxId++;
	}	
}
#else
	int rayId = curMailboxId++;
#endif

	///////////////////////////////////////////////////////////////////////////////
	// Set up 3D DDA for ray
	real NextCrossingT[3], DeltaT[3];
	int Step[3], Out[3], Pos[3];
	
	for (int axis = 0; axis < 3; ++axis) {

		///////////////////////////////////////////////////////////////////////////////
		// Compute current voxel for axis
		Pos[axis] = PosToVoxel(gridIntersect, axis);

		if (r.direction()[axis] >= 0) {
			///////////////////////////////////////////////////////////////////////////////
			// Handle ray with positive direction for voxel stepping
			NextCrossingT[axis] = rayT + (VoxelToPos(Pos[axis]+1, axis) - gridIntersect[axis]) * r.invDirection()[axis];
			DeltaT[axis] = Width[axis] * r.invDirection()[axis];
			Step[axis] = 1;
			Out[axis] = NVoxels[axis];
		}
		else {
			///////////////////////////////////////////////////////////////////////////////
			// Handle ray with negative direction for voxel stepping
			NextCrossingT[axis] = rayT + (VoxelToPos(Pos[axis], axis) - gridIntersect[axis]) * r.invDirection()[axis];
			DeltaT[axis] = -Width[axis] * r.invDirection()[axis];
			Step[axis] = -1;
			Out[axis] = -1;
		}
	}

	///////////////////////////////////////////////////////////////////////////////
	// Walk ray through voxel grid
	bool hitSomething = false;
	while(true) {
		Voxel *voxel = voxels[Offset(Pos[0], Pos[1], Pos[2])];
		if (voxel != NULL)
			hitSomething |= voxel->hit(r, intersection, rayId);

		///////////////////////////////////////////////////////////////////////////////
		// Advance to next voxel
		///////////////////////////////////////////////////////////////////////////////

		///////////////////////////////////////////////////////////////////////////////
		// Find _stepAxis_ for stepping to next voxel
		int bits =	((NextCrossingT[0] < NextCrossingT[1]) << 2) +
					((NextCrossingT[0] < NextCrossingT[2]) << 1) +
					((NextCrossingT[1] < NextCrossingT[2]));

		const int cmpToAxis[8] = { 2, 1, 2, 1, 2, 2, 0, 0 };

		int stepAxis = cmpToAxis[bits];
		if (r.tMax() < NextCrossingT[stepAxis])
			break;

		Pos[stepAxis] += Step[stepAxis];
		if (Pos[stepAxis] == Out[stepAxis])
			break;

		NextCrossingT[stepAxis] += DeltaT[stepAxis];
	}

	return hitSomething;
}

const bool Voxel::hit(const Ray& r, Intersection &intersection, int rayId) {

	///////////////////////////////////////////////////////////////////////////////
	// Refine primitives in voxel if needed
	if (!allCanIntersect) {

		MailboxPrim **mpp;
		if (nPrimitives == 1) mpp = &onePrimitive;
		else mpp = primitives;

		for (unsigned int i = 0; i < nPrimitives; ++i) {
			MailboxPrim *mp = mpp[i];

			///////////////////////////////////////////////////////////////////////////////
			// Refine primitive in mp if it's not intersectable
			if (!mp->primitive.canIntersect()) {
				std::vector<Primitive> p;
				mp->primitive.fullyRefine(p);
				
				if (p.size() == 1)	mp->primitive = p[0];
				else				mp->primitive = GridAccel(p, true, false);
			}
		}
		allCanIntersect = true;
	}

	///////////////////////////////////////////////////////////////////////////////
	// Loop over primitives in voxel and find intersections
	bool hitSomething = false;

	MailboxPrim **mpp;
	if (nPrimitives == 1) mpp = &onePrimitive;
	else mpp = primitives;

	for (unsigned int i = 0; i < nPrimitives; ++i) {
		MailboxPrim *mp = mpp[i];

		///////////////////////////////////////////////////////////////////////////////
		// Do mailbox check between ray and primitive
		if (mp->lastMailboxId == rayId)
			continue;

		///////////////////////////////////////////////////////////////////////////////
		// Check for ray--primitive intersection
		mp->lastMailboxId = rayId;
		
		if (mp->primitive.hit(r, intersection)) {			
			hitSomething = true;
			r.setTmax(intersection.tHit);
		}
	}

	return hitSomething;
}

const bool GridAccel::GridAccelValue::shadowHit(const Ray& r) const {

	///////////////////////////////////////////////////////////////////////////////
	// Check ray against overall grid bounds
	real rayT, tmax;

	if (bounds.isInside(r(r.tMin())))			rayT = r.tMin();
	else if (!bounds.intersectP(r, rayT, tmax))	return false;

	Vector3<real> gridIntersect = r(rayT);

	///////////////////////////////////////////////////////////////////////////////
	// Get ray mailbox id
//	int rayId = ++curMailboxId;
#ifdef USE_OPENMP
	int rayId;
{
#pragma omp critical 
		{
			rayId = curMailboxId++;
		}	
}	
#else
	int rayId = curMailboxId++;
#endif

	///////////////////////////////////////////////////////////////////////////////
	// Set up 3D DDA for ray
	real NextCrossingT[3], DeltaT[3];
	int Step[3], Out[3], Pos[3];

	for (int axis = 0; axis < 3; ++axis) {

		///////////////////////////////////////////////////////////////////////////////
		// Compute current voxel for axis
		Pos[axis] = PosToVoxel(gridIntersect, axis);
		if (r.direction()[axis] >= 0) {
			///////////////////////////////////////////////////////////////////////////////
			// Handle ray with positive direction for voxel stepping
			NextCrossingT[axis] = rayT + (VoxelToPos(Pos[axis]+1, axis) - gridIntersect[axis]) * r.invDirection()[axis];
			DeltaT[axis] = Width[axis] * r.invDirection()[axis];
			Step[axis] = 1;
			Out[axis] = NVoxels[axis];
		}
		else {
			///////////////////////////////////////////////////////////////////////////////
			// Handle ray with negative direction for voxel stepping
			NextCrossingT[axis] = rayT + (VoxelToPos(Pos[axis], axis) - gridIntersect[axis]) * r.invDirection()[axis];
			DeltaT[axis] = -Width[axis] * r.invDirection()[axis];
			Step[axis] = -1;
			Out[axis] = -1;
		}
	}

	///////////////////////////////////////////////////////////////////////////////
	// Walk ray through voxel grid
	while(true) {
		Voxel *voxel = voxels[Offset(Pos[0], Pos[1], Pos[2])];
		if (voxel && voxel->shadowHit(r, rayId))
			return true;

		///////////////////////////////////////////////////////////////////////////////
		// Advance to next voxel

		///////////////////////////////////////////////////////////////////////////////
		// Find _stepAxis_ for stepping to next voxel
		int bits =	((NextCrossingT[0] < NextCrossingT[1]) << 2) +
					((NextCrossingT[0] < NextCrossingT[2]) << 1) +
					((NextCrossingT[1] < NextCrossingT[2]));

		const int cmpToAxis[8] = { 2, 1, 2, 1, 2, 2, 0, 0 };

		int stepAxis = cmpToAxis[bits];
		if (r.tMax() < NextCrossingT[stepAxis])
			break;

		Pos[stepAxis] += Step[stepAxis];
		if (Pos[stepAxis] == Out[stepAxis])
			break;

		NextCrossingT[stepAxis] += DeltaT[stepAxis];
	}

	return false;
}

const bool Voxel::shadowHit(const Ray& r, int rayId) {

	///////////////////////////////////////////////////////////////////////////////
	// Refine primitives in voxel if needed
	if (!allCanIntersect) {

		MailboxPrim **mpp;
		if (nPrimitives == 1) mpp = &onePrimitive;
		else mpp = primitives;

		for (unsigned int i = 0; i < nPrimitives; ++i) {
			MailboxPrim *mp = mpp[i];

			///////////////////////////////////////////////////////////////////////////////
			// Refine primitive in mp if it's not intersectable
			if (!mp->primitive.canIntersect()) {
				std::vector<Primitive> p;
				mp->primitive.fullyRefine(p);

				if (p.size() == 1)	mp->primitive = p[0];
				else				mp->primitive = GridAccel(p, true, false);
			}
		}
		allCanIntersect = true;
	}

	///////////////////////////////////////////////////////////////////////////////
	// Loop over primitives in voxel and find intersections
	MailboxPrim **mpp;
	if (nPrimitives == 1) mpp = &onePrimitive;
	else mpp = primitives;

	for (unsigned int i = 0; i < nPrimitives; ++i) {
		MailboxPrim *mp = mpp[i];

		///////////////////////////////////////////////////////////////////////////////
		// Do mailbox check between ray and primitive
		if (mp->lastMailboxId == rayId)
			continue;

		///////////////////////////////////////////////////////////////////////////////
		// Check for ray_primitive intersection
		mp->lastMailboxId = rayId;

		if (mp->primitive.shadowHit(r)) 
			return true;
	}

	return false;
}
