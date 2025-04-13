#ifndef	_GEOMETRYSPACE_LEVELSET_FMMHEAP_GENERAL_H_
#define	_GEOMETRYSPACE_LEVELSET_FMMHEAP_GENERAL_H_

#include <vector>
#include "../../Arithmetic/GridLibrary/GridLibrary.h"

namespace geometryspace {

template<typename T, typename INDEX = arithmetic::grid_index_type>
class NodeHeap {
public : 
	typedef INDEX index_type;

	index_type		idx;
	T				distance;

public : 
	NodeHeap();
	NodeHeap(index_type, T);
};

template<typename T, typename INDEX, typename CONTAINER>
class FMMHeap_General {
public :
	typedef INDEX index_type;	

	enum { MAX_DISTANCE = 100000 };
	enum { FMM_FAR = 0, FMM_ACCEPTED, FMM_TRIAL	};

public:
	FMMHeap_General(index_type minHeapSize);
	void init(index_type minHeapSize);

	inline const bool empty() const	{ return (heapSize == 0); }	
	const NodeHeap<T,INDEX> pop();

	void push(index_type i, T distance);	
	inline void assignDistance(index_type i, T distance) { Heap_Decrease_Key(nodeToHeapIdx[i], distance);}

protected :
	inline const std::size_t Parent(std::size_t i) const	{ return (i/2);		}
	inline const std::size_t Right(std::size_t i) const		{ return (i*2+1);	}
	inline const std::size_t Left(std::size_t i) const		{ return (i*2);		}		

	inline const NodeHeap<T,INDEX> top() const				{ return (nodeHeapArray[1]); 	}

	inline void Exchange(std::size_t i, std::size_t j);
	void Heap_Decrease_Key(std::size_t i, T newDistance);
	void Min_Heapify(std::size_t i);	

	std::vector<NodeHeap<T,INDEX>>	nodeHeapArray;
	CONTAINER						nodeToHeapIdx;

	std::size_t						heapSize;	
};

template<typename T, typename INDEX> NodeHeap<T,INDEX>::NodeHeap() : idx(-1), distance(FMMHeap_General<T,INDEX>::MAX_DISTANCE)	{}
template<typename T, typename INDEX> NodeHeap<T,INDEX>::NodeHeap(index_type i, T dist) : idx(i), distance(dist)	{}

template<typename T, typename INDEX, typename CONTAINER> 
FMMHeap_General<T,INDEX,CONTAINER>::FMMHeap_General(index_type minHeapSize)	{
	assert(minHeapSize > 0);	

	heapSize = -1;
	init(minHeapSize);
}

template<typename T, typename INDEX, typename CONTAINER> 
void FMMHeap_General<T,INDEX,CONTAINER>::init(index_type minHeapSize) {
	nodeHeapArray.clear();	
	nodeHeapArray.reserve(minHeapSize);		
	nodeHeapArray.push_back(NodeHeap<T,INDEX>(-1, MAX_DISTANCE));
	heapSize = 0;
}

template<typename T, typename INDEX, typename CONTAINER> 
inline void FMMHeap_General<T,INDEX,CONTAINER>::Exchange(std::size_t i, std::size_t j)	{	
	std::swap(nodeHeapArray[i], nodeHeapArray[j]);

	nodeToHeapIdx[nodeHeapArray[i].idx] = i;
	nodeToHeapIdx[nodeHeapArray[j].idx] = j;
}

template<typename T, typename INDEX, typename CONTAINER> 
void FMMHeap_General<T,INDEX,CONTAINER>::Heap_Decrease_Key(std::size_t i, T newDistance)	{
	std::size_t idx = i, j;
	nodeHeapArray[idx].distance = newDistance;
	while(idx > 1 && nodeHeapArray[j = Parent(idx)].distance > nodeHeapArray[idx].distance)	{
		Exchange(idx, j);
		idx = j;
	}
}

template<typename T, typename INDEX, typename CONTAINER> 
void FMMHeap_General<T,INDEX,CONTAINER>::Min_Heapify(std::size_t i) {	
	std::size_t l = Left(i);
	std::size_t r = Right(i);
	std::size_t smallest;

	if(l <= heapSize && nodeHeapArray[l].distance < nodeHeapArray[i].distance)
		smallest = l;
	else 
		smallest = i;

	if(r <= heapSize && nodeHeapArray[r].distance < nodeHeapArray[smallest].distance)
		smallest = r;

	if(smallest != i)	{
		Exchange(i,smallest);
		Min_Heapify(smallest);
	}
}

template<typename T, typename INDEX, typename CONTAINER> 
const NodeHeap<T,INDEX> FMMHeap_General<T,INDEX,CONTAINER>::pop()	{
	NodeHeap<T,INDEX> _min(nodeHeapArray[1]);	

	nodeHeapArray[1] = nodeHeapArray[heapSize--];

	if(heapSize != 0)	
		nodeToHeapIdx[nodeHeapArray[1].idx] = 1;

	if(heapSize > 1)		
		Min_Heapify(1);

	return _min;
}

template<typename T, typename INDEX, typename CONTAINER> 
void FMMHeap_General<T,INDEX,CONTAINER>::push(index_type i, T distance)	{
	const NodeHeap<T,INDEX> node(i, distance);	

	if(++heapSize >= nodeHeapArray.size())	nodeHeapArray.push_back(node);
	else									nodeHeapArray[heapSize] = node;

	nodeToHeapIdx[node.idx] = heapSize;
	Heap_Decrease_Key(heapSize,node.distance);
}

};		// namespace geometryspace

#endif
