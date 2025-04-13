#include <stdafx.h>
#include <Illumination/IlluminationHeader/Primitive/KdTreeAccel.h>

using namespace illuminationspace;

void KdTreeAccel::KdTreeAccelValue::buildTree(int nodeNum, const BBox &nodeBounds, const std::vector<BBox> &allPrimBounds, int *primNums, int nPrims, int depth, BoundEdge *edges[3], int *prims0, int *prims1, int badRefines) {

	///////////////////////////////////////////////////////////////////////////////
	// Get next free node from _nodes_ array
	if (nextFreeNode == nAllocedNodes) {
		int nAlloc = MAX(2 * nAllocedNodes, 512);
		KdAccelNode *n = (KdAccelNode *)basis::AllocAligned(nAlloc * sizeof(KdAccelNode));
		if (nAllocedNodes > 0) {
			memcpy(n, nodes, nAllocedNodes * sizeof(KdAccelNode));
			basis::FreeAligned(nodes);
		}
		nodes = n;
		nAllocedNodes = nAlloc;
	}
	++nextFreeNode;

	///////////////////////////////////////////////////////////////////////////////
	// Initialize leaf node if termination criteria met
	if (nPrims <= maxPrims || depth == 0) {
		nodes[nodeNum].initLeaf(primNums, nPrims, mailboxPrims, arena);
		return;
	}

	///////////////////////////////////////////////////////////////////////////////
	// Initialize interior node and continue recursion
	///////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////
	// Choose split axis position for interior node
	int bestAxis = -1, bestOffset = -1;
	real bestCost = FLT_MAX;
	real oldCost = isectCost * real(nPrims);
	Vector3<real> d = nodeBounds.max() - nodeBounds.min();
	real totalSA = (2.f * (d.x()*d.y() + d.x()*d.z() + d.y()*d.z()));
	real invTotalSA = 1.f / totalSA;

	///////////////////////////////////////////////////////////////////////////////
	// Choose which axis to split along
	int axis;
	if (d.x() > d.y() && d.x() > d.z()) axis = 0;
	else axis = (d.y() > d.z()) ? 1 : 2;

	int retries = 0;
retrySplit:
	///////////////////////////////////////////////////////////////////////////////
	// Initialize edges for _axis_
	for (int i = 0; i < nPrims; ++i) {
		int pn = primNums[i];
		const BBox &bbox = allPrimBounds[pn];
		edges[axis][2*i  ] = BoundEdge(bbox.min()[axis], pn, true);
		edges[axis][2*i+1] = BoundEdge(bbox.max()[axis], pn, false);
	}
	std::sort(&edges[axis][0], &edges[axis][2*nPrims]);

	///////////////////////////////////////////////////////////////////////////////
	// Compute cost of all splits for selected axis to find best
	int nBelow = 0, nAbove = nPrims;
	for (int i = 0; i < 2*nPrims; ++i) {

		if (edges[axis][i].type == BoundEdge::END) --nAbove;
		real edget = edges[axis][i].t;

		if (edget > nodeBounds.min()[axis] && edget < nodeBounds.max()[axis]) {

			///////////////////////////////////////////////////////////////////////////////
			// Compute cost for split at i-th edge
			int otherAxis[3][2] = { {1,2}, {0,2}, {0,1} };
			int otherAxis0 = otherAxis[axis][0];
			int otherAxis1 = otherAxis[axis][1];
			real belowSA = 2 * (d[otherAxis0] * d[otherAxis1] + (edget - nodeBounds.min()[axis]) * (d[otherAxis0] + d[otherAxis1]));
			real aboveSA = 2 * (d[otherAxis0] * d[otherAxis1] + (nodeBounds.max()[axis] - edget) * (d[otherAxis0] + d[otherAxis1]));
			real pBelow = belowSA * invTotalSA;
			real pAbove = aboveSA * invTotalSA;
			real eb = (nAbove == 0 || nBelow == 0) ? emptyBonus : 0.f;
			real cost = traversalCost + isectCost * (1.f - eb) * (pBelow * nBelow + pAbove * nAbove);

			///////////////////////////////////////////////////////////////////////////////
			// Update best split if this is lowest cost so far
			if (cost < bestCost)  {
				bestCost = cost;
				bestAxis = axis;
				bestOffset = i;
			}
		}
		if (edges[axis][i].type == BoundEdge::START) ++nBelow;
	}

	//Assert(nBelow == nPrims && nAbove == 0); // NOBOOK

	///////////////////////////////////////////////////////////////////////////////
	// Create leaf if no good splits were found
	if (bestAxis == -1 && retries < 2) {
		++retries;
		axis = (axis+1) % 3;
		goto retrySplit;
	}
	if (bestCost > oldCost) ++badRefines;
	if ((bestCost > 4.f * oldCost && nPrims < 16) || bestAxis == -1 || badRefines == 3) {
		nodes[nodeNum].initLeaf(primNums, nPrims, mailboxPrims, arena);
		return;
	}

	///////////////////////////////////////////////////////////////////////////////
	// Classify primitives with respect to split
	int n0 = 0, n1 = 0;
	for (int i = 0; i < bestOffset; ++i)
		if (edges[bestAxis][i].type == BoundEdge::START)
			prims0[n0++] = edges[bestAxis][i].primNum;

	for (int i = bestOffset+1; i < 2*nPrims; ++i)
		if (edges[bestAxis][i].type == BoundEdge::END)
			prims1[n1++] = edges[bestAxis][i].primNum;

	///////////////////////////////////////////////////////////////////////////////
	// Recursively initialize children nodes
	real tsplit = edges[bestAxis][bestOffset].t;
	nodes[nodeNum].initInterior(bestAxis, tsplit);

	BBox bounds0 = nodeBounds, bounds1 = nodeBounds;
	bounds0.max()[bestAxis] = bounds1.min()[bestAxis] = tsplit;

	buildTree(nodeNum+1, bounds0, allPrimBounds, prims0, n0, depth-1, edges, prims0, prims1 + nPrims, badRefines);
	nodes[nodeNum].aboveChild = nextFreeNode;
	buildTree(nodes[nodeNum].aboveChild, bounds1, allPrimBounds, prims1, n1, depth-1, edges, prims0, prims1 + nPrims, badRefines);
}