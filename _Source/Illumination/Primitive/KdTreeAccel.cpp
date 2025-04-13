#include <stdafx.h>
#include <Illumination/IlluminationHeader/Illumination/BasicLibrary.h>
#include <Illumination/IlluminationHeader/Primitive/KdTreeAccel.h>

using namespace illuminationspace;

KdTreeAccel::KdTreeAccelValue::KdTreeAccelValue(const std::vector<Primitive> &p, int icost, int tcost, real ebonus, int maxp, int maxDepth)
: isectCost(icost), traversalCost(tcost), maxPrims(maxp), emptyBonus(ebonus) {

	std::vector<Primitive> prims;
	for(unsigned int i=0;i<p.size();++i)
		p[i].fullyRefine(prims);

	///////////////////////////////////////////////////////////////////////////////
	// Initialize MailBox
	curMailboxId = 0;
	nMailboxes = static_cast<unsigned int>(prims.size());
	mailboxPrims = (MailboxPrim*)basis::AllocAligned(nMailboxes * sizeof(MailboxPrim));

	for(unsigned int i=0;i<nMailboxes;++i)
		new (&mailboxPrims[i]) MailboxPrim(prims[i]);

	///////////////////////////////////////////////////////////////////////////////
	// Build KdTree for accelerator
	nextFreeNode = nAllocedNodes = 0;

	if(maxDepth <= 0)
		maxDepth = Round2Int(8 + 1.3f * Log2Int(float(prims.size())));

	///////////////////////////////////////////////////////////////////////////////
	// Compute Bounds for KdTree construction
	std::vector<BBox> primBounds;
	primBounds.reserve(prims.size());
	for(unsigned int i=0;i<prims.size();++i) {
		BBox b = prims[i].worldBound();
		bounds.surround(b);
		primBounds.push_back(b);
	}

	///////////////////////////////////////////////////////////////////////////////
	// Allocate working memory for kdTree construction
	BoundEdge *edges[3];
	for(int i=0;i<3;++i)
		edges[i] = new BoundEdge[2*prims.size()];

	int *prims0 = new int[prims.size()];
	int *prims1 = new int[(maxDepth+1) * prims.size()];

	///////////////////////////////////////////////////////////////////////////////
	// Initialize primNums for kdTree construction
	int *primNums = new int[prims.size()];
	for(unsigned int i=0;i<prims.size();++i)
		primNums[i] = i;

	///////////////////////////////////////////////////////////////////////////////
	// Start recursive construction of KD-tree
	buildTree(0, bounds, primBounds, primNums, static_cast<int>(prims.size()), maxDepth, edges, prims0, prims1);

	///////////////////////////////////////////////////////////////////////////////
	// Free working memory for KD-tree construction	
	for (int i = 0; i < 3; ++i)
		delete[] edges[i];
	delete[] primNums;
	delete[] prims0;
	delete[] prims1;
}