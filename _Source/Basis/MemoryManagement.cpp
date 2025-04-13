#include <stdafx.h>
#include <Basis/MemoryManagement/MemoryManagement.h>

namespace basis {

void *AllocAligned(size_t size) {
#ifndef L1_CACHE_LINE_SIZE
#define L1_CACHE_LINE_SIZE 64
#endif
	return memalign(L1_CACHE_LINE_SIZE, size);
}

void FreeAligned(void *ptr) {
#ifdef WIN32 
	_aligned_free(ptr);
#else 
	free(ptr);
#endif 
}

};	// basis