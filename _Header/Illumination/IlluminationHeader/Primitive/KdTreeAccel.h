#ifndef _ILLUMINATION_KDTREE_ACCEL_H_
#define _ILLUMINATION_KDTREE_ACCEL_H_

#include <limits>
#include <algorithm>
#include "../Illumination/MemoryArena.h"
#include "../Core/BBox.h"
#include "Aggregate.h"

namespace illuminationspace {

struct KdAccelNode;
struct BoundEdge;
class ILLUMINATIONDLL KdTreeAccel : public Aggregate {
protected : 
	class KdTreeAccelValue;

public : 
	KdTreeAccel(const std::vector<Primitive> &p, int icost=80, int tcost=1, real ebonus=0.5f, int maxp=1, int maxDepth=-1)
	: Aggregate(new KdTreeAccelValue(p, icost, tcost, ebonus, maxp, maxDepth)) {}

	void buildTree(int nodeNum, const BBox &nodeBounds, const std::vector<BBox> &allPrimBounds, int *primNums, int nPrims, int depth, BoundEdge *edges[3], int *prims0, int *prims1, int badRefines=0) {
		dynamic_cast<KdTreeAccelValue*>(const_cast<PrimitiveValue*>(value.getptr()))->buildTree(nodeNum, nodeBounds, allPrimBounds, primNums, nPrims, depth, edges, prims0, prims1, badRefines);
	}

protected : 
	class ILLUMINATIONDLL KdTreeAccelValue : public Aggregate::AggregateValue {
	public : 
		KdTreeAccelValue(const std::vector<Primitive> &p, int icost=80, int tcost=1, real ebonus=0.5f, int maxp=1, int maxDepth=-1);

		const BBox worldBound() const { return bounds; }
		const bool canIntersect() const { return true; }
		void buildTree(int nodeNum, const BBox &nodeBounds, const std::vector<BBox> &allPrimBounds, int *primNums, int nPrims, int depth, BoundEdge *edges[3], int *prims0, int *prims1, int badRefines=0);

		const bool hit(const Ray& r, Intersection &intersection) const;
		const bool shadowHit(const Ray& r) const;

	public :
		///////////////////////////////////////////////////////////////////////////////
		// For Construction
		int isectCost, traversalCost, maxPrims;
		real emptyBonus;

		///////////////////////////////////////////////////////////////////////////////
		// For Traversal _ MailBox _ preventing intersection test for same primitives
		unsigned int nMailboxes;
		MailboxPrim *mailboxPrims;
		mutable int curMailboxId;

		///////////////////////////////////////////////////////////////////////////////
		// KD-Tree Node
		KdAccelNode *nodes;
		int nAllocedNodes, nextFreeNode;

		///////////////////////////////////////////////////////////////////////////////
		// Bounding Box
		BBox bounds;

		///////////////////////////////////////////////////////////////////////////////
		// For Memory Management
		MemoryArena arena;
	};
};

struct KdAccelNode {
public :
	void initLeaf(int *primNums, int np, MailboxPrim *mailboxPrims, MemoryArena &arena);
	void initInterior(int axis, real s) {
		split = s;
		flags &= ~3;
		flags |= axis;
	}
	const real splitPos() const { return split; }
	const int nPrimitives() const { return nPrims >> 2; }
	const int splitAxis() const { return flags & 3; }
	const bool isLeaf() const { return (flags & 3) == 3; }

public :
	union {
		unsigned int flags;
		real split;
		unsigned int nPrims;
	};
	union {
		unsigned int aboveChild;
		MailboxPrim *onePrimitive;
		MailboxPrim **primitives;
	};
};

struct BoundEdge {
public :
	BoundEdge() {}
	BoundEdge(real tt, int pn, bool starting) {
		t = tt;
		primNum = pn;
		type = starting ? START : END;
	}
	const bool operator<(const BoundEdge &e) const {
		if (t == e.t)
			return static_cast<int>(type) < static_cast<int>(e.type);
		else return t < e.t;
	}

public :
	real t;
	int primNum;
	enum {START, END} type;
};

};		// illuminationspace

#endif