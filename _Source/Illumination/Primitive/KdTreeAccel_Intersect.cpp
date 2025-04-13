#include <stdafx.h>
#include <Illumination/IlluminationHeader/Primitive/KdTreeAccel.h>

using namespace illuminationspace;

#define MAX_TODO 64
struct KdToDo {
	const KdAccelNode *node;
	real tmin, tmax;
};

const bool KdTreeAccel::KdTreeAccelValue::hit(const Ray& ray, Intersection &intersection) const {

	///////////////////////////////////////////////////////////////////////////////
	// Compute initial parametric range of ray inside kd-tree extent
	real tmin, tmax;
	if (!bounds.intersectP(ray, tmin, tmax))
		return false;

	// Prepare to traverse kdTree for ray
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

	Vector3<real> invDir(ray.invDirection());

	KdToDo todo[MAX_TODO];
	int todoPos = 0;

	///////////////////////////////////////////////////////////////////////////////
	// Traverse KD-tree nodes in order for ray
	bool hit = false;
	const KdAccelNode *node = &nodes[0];
	while (node != NULL) {

		///////////////////////////////////////////////////////////////////////////////
		// Bail out if we found a hit closer than the current node
		if (ray.tMax() < tmin) 
			break;

		if (!node->isLeaf()) {
			///////////////////////////////////////////////////////////////////////////////
			// Process kd-tree interior node
			///////////////////////////////////////////////////////////////////////////////

			///////////////////////////////////////////////////////////////////////////////
			// Compute parametric distance along ray to split plane
			int axis = node->splitAxis();
			real tplane = (node->splitPos() - ray.origin()[axis]) * invDir[axis];

			// Get node children pointers for ray
			const KdAccelNode *firstChild, *secondChild;
			bool belowFirst = (ray.origin()[axis] <= node->splitPos());
			if (belowFirst) {
				firstChild = node + 1;
				secondChild = &nodes[node->aboveChild];
			}
			else {
				firstChild = &nodes[node->aboveChild];
				secondChild = node + 1;
			}

			///////////////////////////////////////////////////////////////////////////////
			// Advance to next child node, possibly enqueue other child
			if (tplane > tmax || tplane < 0)
				node = firstChild;
			else if (tplane < tmin)
				node = secondChild;
			else {
				///////////////////////////////////////////////////////////////////////////////
				// Enqueue _secondChild_ in todo list
				todo[todoPos].node = secondChild;
				todo[todoPos].tmin = tplane;
				todo[todoPos].tmax = tmax;
				++todoPos;
				node = firstChild;
				tmax = tplane;
			}
		}
		else {
			///////////////////////////////////////////////////////////////////////////////
			// Check for intersections inside leaf node
			unsigned int nPrimitives = node->nPrimitives();
			if (nPrimitives == 1) {
				MailboxPrim *mp = node->onePrimitive;

				///////////////////////////////////////////////////////////////////////////////
				// Check one primitive inside leaf node
				if (mp->lastMailboxId != rayId) {
					mp->lastMailboxId = rayId;
					if (mp->primitive.hit(ray, intersection)) {
						hit = true;
						ray.setTmax(intersection.tHit);
					}
				}
			}
			else {
				MailboxPrim **prims = node->primitives;
				for (unsigned int i = 0; i < nPrimitives; ++i) {
					MailboxPrim *mp = prims[i];

					// Check one primitive inside leaf node
					if (mp->lastMailboxId != rayId) {
						mp->lastMailboxId = rayId;
						if (mp->primitive.hit(ray, intersection)) {
							hit = true;
							ray.setTmax(intersection.tHit);
						}
					}
				}
			}

			///////////////////////////////////////////////////////////////////////////////
			// Grab next node to process from todo list
			if (todoPos > 0) {
				--todoPos;
				node = todo[todoPos].node;
				tmin = todo[todoPos].tmin;
				tmax = todo[todoPos].tmax;
			}
			else
				break;
		}
	}

	return hit;
}

const bool KdTreeAccel::KdTreeAccelValue::shadowHit(const Ray& ray) const {

	///////////////////////////////////////////////////////////////////////////////
	// Compute initial parametric range of ray inside kd-tree extent
	real tmin, tmax;
	if (!bounds.intersectP(ray, tmin, tmax))
		return false;

	///////////////////////////////////////////////////////////////////////////////
	// Prepare to traverse kd-tree for rayz
	//int rayId = curMailboxId++;
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

	Vector3<real> invDir(ray.invDirection());

	KdToDo todo[MAX_TODO];
	int todoPos = 0;

	const KdAccelNode *node = &nodes[0];
	while (node != NULL) {

		if (!node->isLeaf()) {

			///////////////////////////////////////////////////////////////////////////////
			// Process kd-tree interior node
			///////////////////////////////////////////////////////////////////////////////

			///////////////////////////////////////////////////////////////////////////////
			// Compute parametric distance along ray to split plane
			int axis = node->splitAxis();
			real tplane = (node->splitPos() - ray.origin()[axis]) * invDir[axis];

			// Get node children pointers for ray
			const KdAccelNode *firstChild, *secondChild;
			bool belowFirst = (ray.origin()[axis] <= node->splitPos());
			if (belowFirst) {
				firstChild = node + 1;
				secondChild = &nodes[node->aboveChild];
			}
			else {
				firstChild = &nodes[node->aboveChild];
				secondChild = node + 1;
			}

			///////////////////////////////////////////////////////////////////////////////
			// Advance to next child node, possibly enqueue other child
			if (tplane > tmax || tplane < 0)
				node = firstChild;
			else if (tplane < tmin)
				node = secondChild;
			else {
				///////////////////////////////////////////////////////////////////////////////
				// Enqueue _secondChild_ in todo list
				todo[todoPos].node = secondChild;
				todo[todoPos].tmin = tplane;
				todo[todoPos].tmax = tmax;
				++todoPos;
				node = firstChild;
				tmax = tplane;
			}			
		}
		else {
			///////////////////////////////////////////////////////////////////////////////
			// Check for shadow ray intersections inside leaf node
			unsigned int nPrimitives = node->nPrimitives();
			if (nPrimitives == 1) {
				MailboxPrim *mp = node->onePrimitive;

				if (mp->lastMailboxId != rayId) {
					mp->lastMailboxId = rayId;
					if (mp->primitive.shadowHit(ray))
						return true;
				}
			}
			else {
				MailboxPrim **prims = node->primitives;
				for (unsigned int i = 0; i < nPrimitives; ++i) {
					MailboxPrim *mp = prims[i];

					if (mp->lastMailboxId != rayId) {
						mp->lastMailboxId = rayId;
						if (mp->primitive.shadowHit(ray))
							return true;
					}
				}
			}
			// Grab next node to process from todo list
			if (todoPos > 0) {
				--todoPos;
				node = todo[todoPos].node;
				tmin = todo[todoPos].tmin;
				tmax = todo[todoPos].tmax;
			}
			else
				break;			
		}
	}

	return false;
}