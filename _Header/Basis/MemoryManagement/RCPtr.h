#ifndef _BASIS_RCPTR_H_
#define _BASIS_RCPTR_H_

#include <iostream>
#include "../Basis.h"
#include "../Basis/Environments.h"

namespace basis {

template<typename T>
class RCPtr {
public : 
	inline RCPtr(T * realPtr = NULL);
	inline RCPtr(const RCPtr& rhs);

	template<typename S> 
	inline RCPtr(const RCPtr<S>& rhs);

	~RCPtr();

	RCPtr& operator=(const RCPtr& rhs);

	template<typename S> 
	RCPtr& operator=(const RCPtr<S>& rhs);

	inline operator bool() const { return pointee != NULL; }

	inline T* operator->() const { return pointee;	}
	inline T& operator*() const	 { return *pointee; }

	inline const T* getptr() const { return pointee; }

private : 
	T* pointee;

	void init();

public :
	friend class RCPtr;
};

template<typename T> 
void RCPtr<T>::init() {
	if(pointee == NULL) return;
	if(pointee->isShareable() == false) {
		//pointee = new T(*pointee);
		pointee = pointee->clone();
	}

	pointee->addReference();
}

template<typename T> inline RCPtr<T>::RCPtr(T * realPtr) : pointee(realPtr) { init(); }

template<typename T> inline RCPtr<T>::RCPtr(const RCPtr<T>& rhs) : pointee(rhs.pointee) { init(); }

template<typename T> template<typename S> 
inline RCPtr<T>::RCPtr(const RCPtr<S>& rhs) : pointee(dynamic_cast<T*>(rhs.pointee)) { init(); }

template<typename T> RCPtr<T>::~RCPtr() { 
	if(pointee) 
		pointee->removeReference();
}

template<typename T> RCPtr<T>& RCPtr<T>::operator=(const RCPtr& rhs) {
	if(pointee != rhs.pointee) {
		T *oldPointee = pointee;
		pointee = rhs.pointee;
		init();
		if(oldPointee) oldPointee->removeReference();
	}
	return *this;
}

template<typename T> template<typename S> RCPtr<T>& RCPtr<T>::operator=(const RCPtr<S>& rhs) {
	if(pointee != rhs.pointee) {
		T *oldPointee = pointee;
		pointee = dynamic_cast<T*>(rhs.pointee);
		init();
		if(oldPointee) oldPointee->removeReference();
	}
	return *this;
}

class BASISDLL RCObject {
public : 
	inline void addReference();
	inline void removeReference();

	inline void markUnshareable();
	inline bool isShareable() const;

	inline bool isShared() const;

	virtual RCObject *clone() const = 0;

	//////////////////////////////////////////////////////////////////////////
	// Test Code
	const int getRefCount() const { return refCount; }

protected :
	inline RCObject();
	inline RCObject(const RCObject& rhs);
	inline RCObject& operator=(const RCObject& rhs);
	virtual ~RCObject() = 0;

private :
	int refCount;
	bool shareable;
};

inline RCObject::RCObject() :refCount(0), shareable(true) {}

inline RCObject::RCObject(const RCObject&) : refCount(0), shareable(true) {}

inline RCObject& RCObject::operator=(const RCObject&) { return *this; }

inline void RCObject::addReference() { 
	++refCount; 
}

inline void RCObject::removeReference() { 
	if(--refCount==0) 
		delete this; 
}

inline void RCObject::markUnshareable() { shareable = false; }

inline bool RCObject::isShareable() const { return shareable; }

inline bool RCObject::isShared() const { return refCount > 1; }

};		// namespace basis

#endif