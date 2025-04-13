#include <stdafx.h>
#include <Illumination/IlluminationHeader/Illumination/MemoryArena.h>
#include <Basis/Basis/BasicLibrary.h>

using namespace basis;
using namespace illuminationspace;

MemoryArena::MemoryArena(unsigned int bs) {
	blockSize = bs;
	curBlockPos = 0;
	currentBlock = (char *)basis::AllocAligned(blockSize);
}

MemoryArena::~MemoryArena() {
	basis::FreeAligned(currentBlock);
	for (unsigned int i = 0; i < usedBlocks.size(); ++i)		basis::FreeAligned(usedBlocks[i]);
	for (unsigned int i = 0; i < availableBlocks.size(); ++i)	basis::FreeAligned(availableBlocks[i]);
}

void *MemoryArena::Alloc(unsigned int sz) {
	// Round up _sz_ to minimum machine alignment
	sz = ((sz + 7) & (~7));
	if (curBlockPos + sz > blockSize) {
		// Get new block of memory for _MemoryArena_
		usedBlocks.push_back(currentBlock);
		if (availableBlocks.size() && sz <= blockSize) {
			currentBlock = availableBlocks.back();
			availableBlocks.pop_back();
		}
		else
			currentBlock = (char *)basis::AllocAligned(MAX(sz, blockSize));
		curBlockPos = 0;
	}
	void *ret = currentBlock + curBlockPos;
	curBlockPos += sz;
	return ret;
}

void MemoryArena::FreeAll() {
	curBlockPos = 0;
	while (usedBlocks.size()) {
		availableBlocks.push_back(usedBlocks.back());
		usedBlocks.pop_back();
	}
}