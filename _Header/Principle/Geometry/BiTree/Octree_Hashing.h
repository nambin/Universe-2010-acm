#pragma once

#include <vector>
#include <stack>
#include "../GeometricPrimitive/AABB.h"

namespace geometryspace {

enum OCTREE_ADD_OPTION { OCTREE_PUSH, OCTREE_REFINE, OCTREE_NONE };

#define _OCTREE_BASE_DEPTH 9

template<typename NodeData, typename T> struct OctNode_Hashing {
public : 
	OctNode_Hashing() {
		for(int i=0;i<8;++i)
			children[i] = NULL;
	}	
	~OctNode_Hashing()	{ release();	}
	void init()			{ release();	}
	void release() {
		for(int i=0;i<8;++i) {
			if(children[i]) 
				delete children[i];
			children[i] = NULL;
		}
		data.clear();
		bound.init();	
	}
	void set(const OctNode_Hashing& node) {
		release();

		for(int i=0;i<8;++i) {
			if(!node.children[i])
				continue;
			children[i] = new OctNode<NodeData>;
			children[i]->Set(*node.children[i]);
		}		
		for(size_t i=0;i<node.data.size();++i)
			data.push_back(node.data[i]);
		bound = node.bound;
	}

public : 
	OctNode_Hashing*		children[8];
	std::vector<NodeData>	data;
	AABB<T>					bound;
};

template <typename NodeData, typename AddupProc, typename T> class Octree_Hashing {
public :
	Octree_Hashing(const AABB<T> &b, int md = _OCTREE_BASE_DEPTH) {
		maxDepth = md;
		root.bound = b;
	}
	~Octree_Hashing()	{ release(); }
	void init()			{ release(); }
	void release()		{ root.init();	}
	void set(const Octree_Hashing *octree) {
		if(octree == NULL)
			return;

		maxDepth = octree->maxDepth;
		root.set(octree->root);
	}
	void add(const NodeData &dataItem, const AddupProc &process) { addPrivate(&root, dataItem, process); }

	template<typename LookupProc> void lookup(const Vector3<T> &p, const LookupProc &process) {	lookupPrivate(&root, p, process); }	

private : 
	int								maxDepth;	
	OctNode_Hashing<NodeData,T>		root;

	void addPrivate(OctNode_Hashing<NodeData,T> *node, const NodeData &dataItem, const AddupProc &process, int depth = 0);
	template<typename LookupProc> void lookupPrivate(OctNode_Hashing<NodeData,T> *node, const Vector3<T> &P, const LookupProc &process);
};

template <class NodeData, typename AddupProc, typename T>
void Octree_Hashing<NodeData, AddupProc, T>::addPrivate( OctNode_Hashing<NodeData, T> *node, const NodeData &dataItem, const AddupProc &process, int depth) {

	const AABB<T> &nodeBound = node->bound;
	OCTREE_ADD_OPTION option = process(dataItem, nodeBound);

	// Possibly add data item to current octree node
	if(option == OCTREE_PUSH || (option == OCTREE_REFINE && depth == maxDepth)) {
		node->data.push_back(dataItem);
		return;
	}
	else if(option == OCTREE_NONE || depth == maxDepth)
		return;

	// Otherwise add data item to octree children
	Vector3<T> pMid = T(0.5f) * (nodeBound.min() + nodeBound.max());

	for (int child = 0; child < 8; ++child) {
		if (!node->children[child])	{
			node->children[child] = new OctNode_Hashing<NodeData, T>;

			// Compute childBound for octree child
			AABB<T> &childBound = node->children[child]->bound;
			childBound.min().x() = (child & 4) ? pMid.x() : nodeBound.min().x();
			childBound.max().x() = (child & 4) ? nodeBound.max().x() : pMid.x();
			childBound.min().y() = (child & 2) ? pMid.y() : nodeBound.min().y();
			childBound.max().y() = (child & 2) ? nodeBound.max().y() : pMid.y();
			childBound.min().z() = (child & 1) ? pMid.z() : nodeBound.min().z();
			childBound.max().z() = (child & 1) ? nodeBound.max().z() : pMid.z();			
		}

		addPrivate(node->children[child], dataItem, process, depth+1);
	}
}

template <class NodeData, typename AddupProc, typename T> template<typename LookupProc> 
void Octree_Hashing<NodeData, AddupProc,T>::lookupPrivate( OctNode_Hashing<NodeData, T> *node, const Vector3<T> &p, const LookupProc &process) {

	const AABB<T>& nodeBound = node->bound;

	// Determine whether this node will be looked
	if(!process(p, nodeBound))
		return;

	for (unsigned int i = 0; i < node->data.size(); ++i)
		process(p, node->data[i]);

	// Determine which octree child node p is inside
	Vector3<T> pMid = T(0.5f) * (nodeBound.min() + nodeBound.max());
	unsigned int child = (p.x() > pMid.x() ? 4 : 0) + (p.y() > pMid.y() ? 2 : 0) + (p.z() > pMid.z() ? 1 : 0);

	if (node->children[child]) 
		lookupPrivate(node->children[child], p, process);

	for(unsigned int i=0;i<8;++i) {
		if(child == i)
			continue;

		if (node->children[i]) 
			lookupPrivate(node->children[i], p, process);
	}

//	std::stack<OctNode_Hashing<NodeData, T>*>	node_ptr;
//	std::stack<int>								processing_state;
//
//	node_ptr.push(_node);
//	processing_state.push(0);
//
//	do 
//	{
//state0 :
//		OctNode_Hashing<NodeData, T> *	node = node_ptr.top();				node_ptr.pop();
//		int								state = processing_state.top();		processing_state.pop();
//
//		const AABB<T>& nodeBound = node->bound;
//
//		// Determine whether this node will be looked
//		if(state == 0 && !process(p, nodeBound))
//			continue;
//		state = 1;
//
//		if(state == 1) {
//			for (unsigned int i = 0; i < node->data.size(); ++i)
//				process(p, node->data[i]);
//		}
//		state = 2;
//
//		// Determine which octree child node p is inside
//		Vector3<T> pMid = T(0.5f) * (nodeBound.min() + nodeBound.max());
//		int child = (p.x() > pMid.x() ? 4 : 0) + (p.y() > pMid.y() ? 2 : 0) + (p.z() > pMid.z() ? 1 : 0);
//
//		if(state == 2 && node->children[child]) {
//			node_ptr.push(node);
//			processing_state.push(state);
//
//			node_ptr.push(node->children[child]);
//			processing_state.push(0);
//
//			goto state0;
//		}
//		state = 3;
//
////state3 :
//		for(int i=0;i<8;++i) {
//			if(child == i)
//				continue;
//
//			if (state-3 <= i && node->children[i]) {
//				node_ptr.push(node);
//				processing_state.push(state);
//
//				node_ptr.push(node->children[child]);
//				processing_state.push(0);
//
//				goto state0;
//			}
//			state = 4+i;
//		}
//
//	} while(node_ptr.empty() == false);

	
}

#undef _OCTREE_BASE_DEPTH

};