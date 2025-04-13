#ifndef _ILLUMINATION_MEMORYARENA_H_
#define _ILLUMINATION_MEMORYARENA_H_

#include <vector>
#include "../../../Basis/MemoryManagement/MemoryManagement.h"
#include "../../IlluminationSpace.h"

namespace illuminationspace {

class ILLUMINATIONDLL MemoryArena {
public:
	///////////////////////////////////////////////////////////////////////////////
	// MemoryArena Public Methods
	MemoryArena(unsigned int bs = 32768);
	~MemoryArena();

	void *Alloc(unsigned int sz);
	void FreeAll();

private :
	unsigned int			curBlockPos, blockSize;
	char*					currentBlock;
	std::vector<char*>		usedBlocks, availableBlocks;
};

};		// illuminationspace

#endif