#ifndef _BASIS_ARRAY_H_
#define _BASIS_ARRAY_H_

#include <cassert>
#include <iostream>
#include <vector>

namespace basis {

#ifdef _DEBUG
#define __ASSERT(c) assert(c)
#else
#define __ASSERT(c) {}
#endif

template <typename T> class Array {
public : 
	typedef std::size_t			index_type ;
	typedef index_type*			index_pointer ;
	typedef const index_type*	index_const_pointer ;
	typedef index_type&			index_reference ;
	typedef const index_type&	index_const_reference ;

	typedef T					value_type ;
	typedef value_type*			value_pointer ;       
	typedef const value_type*	value_const_pointer ; 
	typedef value_type&			value_reference ;     
	typedef const value_type&	value_const_reference;

	//////////////////////////////////////////////////////////////////////////
	// Constructors/Destructors
	Array();
	Array(const Array &src);
	explicit Array(index_type size);
	virtual ~Array();

	//////////////////////////////////////////////////////////////////////////
	// Memory 
	void clear();
	void resize(index_type _size = 0);
	void release();

	//////////////////////////////////////////////////////////////////////////
	// Index operators
	inline const T&	operator[](index_type i) const;
	inline T&		operator[](index_type i);
	inline const T&	operator()(index_type i) const;
	inline T&		operator()(index_type i);
	inline const T&	get(index_type i) const;

	//////////////////////////////////////////////////////////////////////////
	// Assignment and set
	void resize_set(const Array &src);
	void resize_set(const std::vector<T> &src);

	void resize_memcpy(const Array& src);

	void set(const Array &src);
	void set(const T *src);
	void set(const std::vector<T> &src);
	void set(const T &val);

	void memcpy(const Array &src);
	void memsetZero();

	Array &operator=(const Array &src);	
	Array &operator=(const std::vector<T> &src);
	
	//////////////////////////////////////////////////////////////////////////
	// Size
	inline const index_type size() const	{ return _size;		}
	inline const bool empty() const			{ return _size == 0;}

	//////////////////////////////////////////////////////////////////////////
	// In place arithmetic
	Array<T> &operator+=(const Array &_array);
	Array<T> &operator-=(const Array &_array);

	template<typename S> Array<T> &operator*=(S scale);
	template<typename S> Array<T> &operator/=(S scale);

	Array<T>& negate();

	//////////////////////////////////////////////////////////////////////////
	// Manipulation
	Array<T>& sort(index_type iMin=0, index_type iMax=INT_MAX);

	//////////////////////////////////////////////////////////////////////////
	// Get Pointer
	const T* const getptr() const { return data; }
	const bool isNull() const { return data == NULL; }
	operator const bool() const { return data != NULL; }

	//////////////////////////////////////////////////////////////////////////
	// IO
	template <typename T> friend std::ostream& operator<<(std::ostream &os, const Array<T> &src);	

	//////////////////////////////////////////////////////////////////////////
	// Friend Methods
	template <typename T> friend void swap(Array<T> &a, Array<T> &b);		

private : 
	index_type	_size;
	T*			data;
};

//////////////////////////////////////////////////////////////////////////
// Constructors
template <typename T> Array<T>::Array() : _size(0), data(NULL) {}
template <typename T> Array<T>::Array(const Array &src) : _size(src._size), data(NULL) {
	resize(src.size());
	set(src);
}
template <typename T> Array<T>::Array(index_type size) : _size(size), data(NULL) { 
	if(_size != 0)
		data = new T[_size];
}

template <typename T> Array<T>::~Array() {
	release(); 
}

//////////////////////////////////////////////////////////////////////////
// Memory 
template <typename T> void Array<T>::clear() {
	resize(0);
}

template <typename T> void Array<T>::resize(index_type __size) {
	release();		
	_size = __size;

	if(_size != 0) 
		data = new T[_size];
}

template <typename T> void Array<T>::release()	{
	if(data != NULL) 
		delete[] data; 
	data = NULL;
	_size = 0;
}

//////////////////////////////////////////////////////////////////////////
// Index operators
template <typename T> inline const T &Array<T>::operator[](index_type i) const { 
	__ASSERT(i<_size); 
	return data[i];	
}

template <typename T> inline T &Array<T>::operator[](index_type i)	{ 
	__ASSERT(i<_size); 
	return data[i];	
}

template <typename T> inline const T &Array<T>::operator()(index_type i) const { 
	__ASSERT(i<_size); 
	return data[i];	
}

template <typename T> inline T &Array<T>::operator()(index_type i)	{ 
	__ASSERT(i<_size); 
	return data[i];	
}

template <typename T> inline const T &Array<T>::get(index_type i) const	{ 
	__ASSERT(i<_size); 
	return data[i];	
}

//////////////////////////////////////////////////////////////////////////
// Assignment and set
template <typename T> void Array<T>::resize_set(const Array &src)	{ 
	resize(src.size());
	set(src);
}

template <typename T> void Array<T>::resize_set(const std::vector<T> &src)	{ 
	resize(src.size());
	set(src);
}

template <typename T> void Array<T>::resize_memcpy(const Array &src)	{ 
	resize(src.size());
	memcpy(src);
}

template <typename T> void Array<T>::set(const Array &src)	{ 
	__ASSERT(_size == src._size);
	for(index_type i=0;i<_size;++i)
		data[i] = src.data[i];	
}

template <typename T> void Array<T>::set(const T *src)	{ 
	for(index_type i=0;i<_size;++i)
		data[i] = src[i];
}

template <typename T> void Array<T>::set(const std::vector<T> &src)	{ 
	__ASSERT(_size == src.size());
	for(index_type i=0;i<_size;++i)
		data[i] = src[i];	
}

template <typename T> void Array<T>::set(const T &val)	{ 
	for(index_type i=0;i<_size;++i)
		data[i] = val;	
}	

template <typename T> void Array<T>::memcpy(const Array &src)	{ 
	__ASSERT(_size == src._size);
	std::memcpy(data, src.data, _size * sizeof(T));	
}

template <typename T> void Array<T>::memsetZero() {
	std::memset(data, 0, _size * sizeof(T));
}

template <typename T> Array<T> &Array<T>::operator=(const Array &src) { 
	set(src);	
	return (*this);
}

template <typename T> Array<T> &Array<T>::operator=(const std::vector<T> &src) { 
	set(src);	
	return (*this);
}

//////////////////////////////////////////////////////////////////////////
// In place arithmetic
template <typename T> Array<T> &Array<T>::operator+=(const Array &_array) { 
	__ASSERT(_size == _array._size); 
	for(index_type i=0;i<_size;++i) 
		data[i] += _array.data[i];
	return (*this);
}

template <typename T> Array<T> &Array<T>::operator-=(const Array &_array) { 
	__ASSERT(_size == _array._size); 
	for(index_type i=0;i<_size;++i) 
		data[i] -= _array.data[i];
	return (*this);
}

template<typename T> template<typename S> Array<T> &Array<T>::operator*=(S scale) { 
	for(index_type i=0;i<_size;++i)
		data[i] *= scale;
	return (*this);
}	

template<typename T> template<typename S> Array<T> &Array<T>::operator/=(S scale) { 
	const S s = T(1.0) / scale;
	for(index_type i=0;i<_size;++i)
		data[i] *= s;
	return (*this);
}

template <typename T> Array<T> &Array<T>::negate() { 
	return (*this) *= T(-1.0);
}

//////////////////////////////////////////////////////////////////////////
// Manipulation
template <typename T> Array<T> &Array<T>::sort(index_type iMin, index_type iMax) {
	__ASSERT(iMin>=0 && iMax>=iMin);

	if(iMax > _size)
		iMax = _size;
	std::sort(&data[iMin], &data[iMax]);

	return (*this);
}

//////////////////////////////////////////////////////////////////////////
// IO
template <typename T> std::ostream& operator<<(std::ostream &os, const Array<T> &src) {	
	for(index_type i=0;i<src.size();++i) 
		os << src[i] << std::endl;
	return os;
}

//////////////////////////////////////////////////////////////////////////
// SWAP
template <typename T> inline void swap(Array<T> &a, Array<T> &b)	{
	T *tmp = a.data;	a.data = b.data;	b.data = tmp;
}

#undef __ASSERT

};	// namespace basis

#endif

