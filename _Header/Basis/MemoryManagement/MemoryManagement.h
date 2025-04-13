#ifndef _BASIS_MEMORYMANAGEMENT_H_
#define _BASIS_MEMORYMANAGEMENT_H_

#include <malloc.h>
#include <vector>
#include "../Basis.h"

namespace basis {

#if defined(WIN32)
#define memalign(a,b) _aligned_malloc(b, a)
#elif defined(__APPLE__)
#define memalign(a,b) valloc(b)
#endif

BASISDLL void* AllocAligned(std::size_t size);
BASISDLL void  FreeAligned(void *ptr);

};		// namespace basis

#endif