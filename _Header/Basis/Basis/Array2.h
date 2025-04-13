#ifndef _BASIS_ARRAY2_H_
#define _BASIS_ARRAY2_H_

#include <cassert>
#include <iostream>
#include <vector>

namespace basis {

#ifdef _DEBUG
#define __ASSERT(c) assert(c)
#else
#define __ASSERT(c) {}
#endif

template <typename T> class Array2 {
public : 
	typedef std::size_t			index_type;
	typedef index_type*			index_pointer ;
	typedef const index_type*	index_const_pointer ;
	typedef index_type&			index_reference ;
	typedef const index_type&	index_const_reference ;

	typedef T					value_type;
	typedef value_type*			value_pointer ;       
	typedef const value_type*	value_const_pointer ; 
	typedef value_type&			value_reference ;     
	typedef const value_type&	value_const_reference;

	//////////////////////////////////////////////////////////////////////////
	// Constructors/Destructors
	Array2();
	Array2(const Array2 &src);
	Array2(index_type nx, index_type ny);
	~Array2();

	//////////////////////////////////////////////////////////////////////////
	// Memory 
	void clear();
	void resize(index_type nx=0, index_type ny=0);
	void release();

	//////////////////////////////////////////////////////////////////////////
	// Index operators
	inline const index_type idx(index_type i, index_type j) const;

	inline const T&	operator[](index_type i) const;
	inline T&		operator[](index_type i);
	inline const T&	operator()(index_type i) const;
	inline T&		operator()(index_type i);
	inline const T&	operator()(index_type i, index_type j) const;
	inline T&		operator()(index_type i, index_type j);
	inline const T&	get(index_type i, index_type j) const;

	//////////////////////////////////////////////////////////////////////////
	// Assignment and set
	void resize_set(const Array2 &src);
	void resize_memcpy(const Array2 &src);

	void set(const Array2 &src);
	void set(const T *src);
	void set(const std::vector<T> &src);
	void set(const T &val);

	void memcpy(const Array2 &src);

	Array2 &operator=(const Array2 &src);	
	Array2 &operator=(const std::vector<T> &src);

	//////////////////////////////////////////////////////////////////////////
	// Size
	inline const index_type size() const	{ return _size;		}
	inline const index_type dimX() const	{ return _nX;		}
	inline const index_type dimY() const	{ return _nY;		}
	inline const bool empty() const			{ return _size == 0;}

	//////////////////////////////////////////////////////////////////////////
	// In place arithmetic
	Array2<T> &operator+=(const Array2 &_array);
	Array2<T> &operator-=(const Array2 &_array);

	template<typename S> Array2<T> &operator*=(S scale);
	template<typename S> Array2<T> &operator/=(S scale);

	Array2<T> &negate();

	//////////////////////////////////////////////////////////////////////////
	// Get Pointer
	const T* const getptr() const { return data; }
	const bool isNull() const { return data == NULL; }
	operator const bool() const { return data != NULL; }

	//////////////////////////////////////////////////////////////////////////
	// IO
	template <typename T> friend std::ostream& operator<<(std::ostream &os, const Array2<T> &src);	

	//////////////////////////////////////////////////////////////////////////
	// Friend Methods
	template <typename T> friend void swap(Array2<T> &a, Array2<T> &b);		

private : 
	index_type _nX, _nY;
	index_type _size;
	T *data;
};

//////////////////////////////////////////////////////////////////////////
// Constructors
template <typename T> Array2<T>::Array2() : _nX(0), _nY(0), _size(0), data(NULL) {}
template <typename T> Array2<T>::Array2(const Array2 &src) : _nX(0), _nY(0), _size(0), data(NULL) {
	resize(src._nX, src._nY);
	set(src);
}
template <typename T> Array2<T>::Array2(index_type nx, index_type ny) : _nX(0), _nY(0), _size(0), data(NULL) { 
	if(nx != 0 && ny != 0) {
		_nX = nx; _nY = ny; _size = nx*ny;
		data = new T[_size];
	}
}

template <typename T> Array2<T>::~Array2() {
	release(); 
}

//////////////////////////////////////////////////////////////////////////
// Memory 
template <typename T> void Array2<T>::clear() {
	resize(0,0);
}

template <typename T> void Array2<T>::resize(index_type nx, index_type ny) {
	release();		
	
	if(nx != 0 && ny != 0) {
		_nX = nx; _nY = ny; _size = nx*ny;
		data = new T[_size];
	}	
	else if(nx != 0) {
		std::cout << "Error in Array2::resize" << std::endl;
		assert(false);
	}	
}

template <typename T> void Array2<T>::release()	{
	if(data != NULL) 
		delete[] data; 
	data = NULL;
	_nX = _nY = _size = 0;
}

//////////////////////////////////////////////////////////////////////////
// Index operators
template <typename T> inline const T &Array2<T>::operator[](index_type i) const { 
	__ASSERT(i<_size); 
	return data[i];	
}

template <typename T> inline T &Array2<T>::operator[](index_type i)	{ 
	__ASSERT(i<_size); 
	return data[i];	
}

template <typename T> inline const T &Array2<T>::operator()(index_type i) const { 
	__ASSERT(i<_size); 
	return data[i];	
}

template <typename T> inline T &Array2<T>::operator()(index_type i)	{ 
	__ASSERT(i<_size); 
	return data[i];	
}

template <typename T> inline const typename Array2<T>::index_type Array2<T>::idx(index_type i, index_type j) const { 
	__ASSERT(i<_nX && j<_nY); 
	return i*_nY + j;	
}

template <typename T> inline const T &Array2<T>::operator()(index_type i, index_type j) const { 
	__ASSERT(i<_nX && j<_nY); 
	return data[i*_nY + j];	
}

template <typename T> inline T &Array2<T>::operator()(index_type i, index_type j)	{ 
	__ASSERT(i<_nX && j<_nY); 
	return data[i*_nY + j];	
}

template <typename T> inline const T &Array2<T>::get(index_type i, index_type j) const	{ 
	__ASSERT(i<_nX && j<_nY); 
	return data[i*_nY + j];	
}

//////////////////////////////////////////////////////////////////////////
// Assignment and set
template <typename T> void Array2<T>::resize_set(const Array2 &src) {
	resize(src.dimX(), src.dimY());
	set(src);
}

template <typename T> void Array2<T>::resize_memcpy(const Array2 &src) {
	resize(src.dimX(), src.dimY());
	memcpy(src);
}

template <typename T> void Array2<T>::set(const Array2 &src)	{ 
	__ASSERT(_size == src._size);
	for(index_type i=0;i<_size;++i)
		data[i] = src.data[i];
}

template <typename T> void Array2<T>::set(const T *src)	{ 
	for(index_type i=0;i<_size;++i)
		data[i] = src[i];
}

template <typename T> void Array2<T>::set(const std::vector<T> &src)	{ 
	__ASSERT(_size == src.size());
	for(index_type i=0;i<_size;++i)
		data[i] = src[i];
}

template <typename T> void Array2<T>::set(const T &val)	{ 
	for(index_type i=0;i<_size;++i)
		data[i] = val;
}	

template <typename T> void Array2<T>::memcpy(const Array2 &src)	{ 
	__ASSERT(_size == src._size);
	std::memcpy(data, src.data, _size * sizeof(T));
}

template <typename T> Array2<T> &Array2<T>::operator=(const Array2 &src) { 
	set(src);	
	return (*this);
}

template <typename T> Array2<T> &Array2<T>::operator=(const std::vector<T> &src) { 
	set(src);	
	return (*this);
}

//////////////////////////////////////////////////////////////////////////
// In place arithmetic
template <typename T> Array2<T> &Array2<T>::operator+=(const Array2 &_array) { 
	__ASSERT(_size == _array._size); 
	for(index_type i=0;i<_size;++i) 
		data[i] += _array.data[i];
	return (*this);
}

template <typename T> Array2<T> &Array2<T>::operator-=(const Array2 &_array) { 
	__ASSERT(_size == _array._size); 
	for(index_type i=0;i<_size;++i) 
		data[i] -= _array.data[i];
	return (*this);
}

template<typename T> template<typename S> Array2<T> &Array2<T>::operator*=(S scale) { 
	for(index_type i=0;i<_size;++i)
		data[i] *= scale;
	return (*this);
}	

template<typename T> template<typename S> Array2<T> &Array2<T>::operator/=(S scale) { 
	const S s = T(1.0) / scale;
	for(index_type i=0;i<_size;++i)
		data[i] *= s;
	return (*this);
}

template <typename T> Array2<T> &Array2<T>::negate() { 
	return (*this) *= T(-1.0);
}

//////////////////////////////////////////////////////////////////////////
// IO
template <typename T> std::ostream& operator<<(std::ostream &os, const Array2<T> &src) {	
	for(index_type i=0;i<src._nX;++i) {
		for(index_type j=0;j<src._nY;++j) {
			os << src[src.idx(i,j)] << " ";
		}
		os << std::endl;
	}
	return os;
}

//////////////////////////////////////////////////////////////////////////
// SWAP
template <typename T> inline void swap(Array2<T> &a, Array2<T> &b)	{
	T *tmp = a.data;	a.data = b.data;	b.data = tmp;
}

#undef __ASSERT

};		// namespace basis

#endif