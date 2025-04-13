#pragma once

#include <hash_map>
#include "FMMHeap_General.h"

namespace geometryspace {

template<typename T, typename INDEX = arithmetic::grid_index_type, typename CONTAINER = std::vector<std::size_t>>
class FMMHeap : public FMMHeap_General<T,INDEX,CONTAINER> {
public : 
	FMMHeap(index_type minHeapSize, index_type numNodes) : FMMHeap_General(minHeapSize) { nodeToHeapIdx.resize(numNodes); }
};

template<typename T, typename INDEX, typename CONTAINER = typename stdext::hash_map<INDEX, std::size_t>>
class FMMHeap_Hash : public FMMHeap_General<T,INDEX,CONTAINER> {
public : 
	FMMHeap_Hash(index_type minHeapSize) : FMMHeap_General(minHeapSize) {}
};

};			// geometryspace