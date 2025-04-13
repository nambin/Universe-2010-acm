#ifndef _ILLUMINATION_OBJECTARENA_H_
#define _ILLUMINATION_OBJECTARENA_H_

#include <vector>
#include "../../../Basis/MemoryManagement/MemoryManagement.h"

namespace illuminationspace {

template <class T> class ObjectArena {
public:
	// ObjectArena Public Methods
	ObjectArena()	{ nAvailable = 0;	}
	~ObjectArena()	{ FreeAll();		}

	T* Alloc();
	operator T*() {	return Alloc();	}	

	void FreeAll();

private :	
	T*					mem;
	int					nAvailable;
	std::vector<T*>		toDelete;
};

template<typename T> T* ObjectArena<T>::Alloc() {
	if (nAvailable == 0) {
		int nAlloc = MAX((unsigned long)16, (unsigned long)(65536/sizeof(T)));
		mem = (T *)basis::AllocAligned(nAlloc * sizeof(T));
		nAvailable = nAlloc;
		toDelete.push_back(mem);
	}
	--nAvailable;
	return mem++;
}

template<typename T> void ObjectArena<T>::FreeAll() {
	for (unsigned int i = 0; i < toDelete.size(); ++i)
		basis::FreeAligned(toDelete[i]);
	toDelete.erase(toDelete.begin(), toDelete.end());
	nAvailable = 0;
}

};		// illuminationspace

#endif