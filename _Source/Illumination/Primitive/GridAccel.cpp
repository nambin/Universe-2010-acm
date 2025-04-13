#include <stdafx.h>
#include <Basis/MemoryManagement/MemoryManagement.h>
#include <Illumination/IlluminationHeader/Primitive/GridAccel.h>

using namespace illuminationspace;

GridAccel::GridAccelValue::GridAccelValue(const std::vector<Primitive> &p, bool forRefined, bool refineImmediately) 
: gridForRefined(forRefined) {

	std::vector<Primitive> prims;
	if(refineImmediately)
		for(unsigned int i=0;i<p.size();++i)
			p[i].fullyRefine(prims);
	else
		prims = p;

	///////////////////////////////////////////////////////////////////////////////
	// Initialize MailBoxes for Grid
	nMailboxes = static_cast<unsigned int>(prims.size());
	mailboxes = (MailboxPrim *)basis::AllocAligned(nMailboxes * sizeof(MailboxPrim));

	for(unsigned int i=0;i<nMailboxes;++i)
		new (&mailboxes[i]) MailboxPrim(prims[i]);

	///////////////////////////////////////////////////////////////////////////////
	// Compute Bounds and Choose Grid Resolution
	for(unsigned int i=0;i<prims.size();++i)
		bounds.surround(prims[i].worldBound());

	Vector3<real> delta = bounds.max() - bounds.min();

	int maxAxis = bounds.maximumExtent();
	real invMaxWidth = 1.0f / delta[maxAxis];
	real cubeRoot = 3.0f * pow(static_cast<real>(prims.size()), static_cast<real>(1.0f/3.0f));
	real voxelsPerUnitDist = cubeRoot * invMaxWidth;

	for(int axis = 0;axis < 3; ++axis) {
		NVoxels[axis] = Round2Int(delta[axis] * voxelsPerUnitDist);
		NVoxels[axis] = clamp(NVoxels[axis], 1, 64);
	}

	///////////////////////////////////////////////////////////////////////////////
	// Compute Voxel Width and Allocate voxels
	for(int axis=0;axis<3;++axis) {
		Width[axis] = delta[axis] / NVoxels[axis];
		InvWidth[axis] = (Width[axis] == 0.0f) ? 0.0f : 1.0f / Width[axis];		
	}
	int nVoxels = NVoxels[0] * NVoxels[1] * NVoxels[2];
	voxels = (Voxel**)basis::AllocAligned(nVoxels * sizeof(Voxel*));
	memset(voxels, 0, nVoxels * sizeof(Voxel*));

	///////////////////////////////////////////////////////////////////////////////
	// Add primitives to grid voxels
	for(unsigned int i=0;i<prims.size();++i) {
		BBox pb = prims[i].worldBound();
		int vmin[3], vmax[3];
		for(int axis=0;axis<3;++axis) {
			vmin[axis] = PosToVoxel(pb.min(), axis);
			vmax[axis] = PosToVoxel(pb.max(), axis);
		}

		for(int x=vmin[0];x<=vmax[0];++x) 
			for(int y=vmin[1];y<=vmax[1];++y) 
				for(int z=vmin[2];z<=vmax[2];++z) {
					int offset = Offset(x,y,z);
					if(!voxels[offset]) {
						//voxels[offset] = new (voxelArena) Voxel(&mailboxes[i]);
						voxels[offset] = new Voxel(&mailboxes[i]);
					}
					else {
						voxels[offset]->AddPrimitive(&mailboxes[i]);
					}
				}
	}

	curMailboxId = 0;
}

void Voxel::AddPrimitive(MailboxPrim *prim) {

	if(nPrimitives == 1) {
		MailboxPrim **p = new MailboxPrim*[2];
		p[0] = onePrimitive;
		primitives = p;
	}
	else if(IsPowerOf2(nPrimitives)) {
		int nAlloc = 2 * nPrimitives;
		MailboxPrim **p = new MailboxPrim *[nAlloc];
		for (unsigned int i=0;i<nPrimitives;++i) 
			p[i] = primitives[i];
		delete[] primitives;
		primitives = p;
	}

	primitives[nPrimitives] = prim;
	++nPrimitives;
}