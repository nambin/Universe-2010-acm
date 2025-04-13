#ifndef _ARITHMETIC_GRID3x4_H_
#define _ARITHMETIC_GRID3x4_H_

#include "../../../Basis/Basis/Environments.h"
#include "../../../Basis/Basis/BasicLibrary.h"
#include "../../../Basis/ParallelProgramming/ParallelComputable.h"
#include "../../../Basis/MemoryManagement/MemoryManagement.h"
#include "GridLibrary.h"
#include "Grid3.h"

namespace arithmetic {

#ifdef _DEBUG
#define __ASSERT(c) assert(c)
#else
#define __ASSERT(c) {}
#endif

template <typename T> class Grid3x4	{
public : 
	typedef typename Grid3<T,T>::index_type	index_type;
	typedef index_type*						index_pointer ;
	typedef const index_type*				index_const_pointer ;
	typedef index_type&						index_reference ;
	typedef const index_type&				index_const_reference ;

	typedef T					value_type;
	typedef value_type*			value_pointer ;       
	typedef const value_type*	value_const_pointer ; 
	typedef value_type&			value_reference ;     
	typedef const value_type&	value_const_reference;

	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Grid3x4();
	Grid3x4(index_type nX, index_type nY, index_type nZ);
	~Grid3x4();

	void init(index_type nX, index_type nY, index_type nZ);
	void set(index_type nX, index_type nY, index_type nZ);
	void release();

	index_type size() const { return _size; }

	//////////////////////////////////////////////////////////////////////////
	// Index operators
	inline const T	operator[](index_type i) const;
	inline T&		operator[](index_type i)		;
	inline const T	operator()(index_type i) const;
	inline T&		operator()(index_type i)		;
	inline const T	operator()(index_type i, index_type j, index_type k, index_type l) const;
	inline T&		operator()(index_type i, index_type j, index_type k, index_type l)		;

	inline const T	get(index_type i, index_type j, index_type k, index_type l) const;

	//////////////////////////////////////////////////////////////////////////
	// Get Pointer
	T* getptr() const { return data; }

	//////////////////////////////////////////////////////////////////////////
	// Assignment and set
	void zero();

	void set(const Grid3x4 &grid);
	void set(const T &val);

	const Grid3x4 &operator=(const Grid3x4 &grid);

	//////////////////////////////////////////////////////////////////////////
	// Manipulation	
	const T infNorm() const;
	const T l2Norm() const;

	const Grid3x4 &operator*=(const T c);

public : 
	index_type	nX, nY, nZ;
	index_type	_size;
	T *			data;
};

//////////////////////////////////////////////////////////////////////////
// Constructor/Destructor
template<typename T> Grid3x4<T>::Grid3x4()
: nX(0), nY(0), nZ(0), _size(0), data(NULL)
{}

template<typename T> Grid3x4<T>::Grid3x4(index_type _nX, index_type _nY, index_type _nZ)
: nX(0), nY(0), nZ(0), _size(0), data(NULL) {
	set(_nX, _nY, _nZ); 
}

template<typename T> Grid3x4<T>::~Grid3x4() { 
	release(); 
}

template<typename T> void Grid3x4<T>::init(index_type _nX, index_type _nY, index_type _nZ)	{
	set(_nX, _nY, _nZ);
}

template<typename T> void Grid3x4<T>::set(index_type _nX, index_type _nY, index_type _nZ)	{
	nX = _nX; nY = _nY; nZ = _nZ;

	release();
	nX = _nX;
	nY = _nY;
	nZ = _nZ;
	_size = nX * nY * nZ * 4;

	data = new T[_size];
	//data = (T*)memoryspace::AllocAligned(_size*sizeof(T));

	zero();
}

template<typename T> void Grid3x4<T>::release()	{
	if(data) delete[] data; data = NULL;
	//if(data) memoryspace::FreeAligned(data); data = NULL;
	nX = nY = nZ = _size = 0;
}

//////////////////////////////////////////////////////////////////////////
// Index operators
template<typename T> inline const T Grid3x4<T>::operator[](index_type i) const	{ __ASSERT(i < _size); return data[i]; }
template<typename T> inline const T Grid3x4<T>::operator()(index_type i) const  { __ASSERT(i < _size); return data[i]; }

template<typename T> inline T &Grid3x4<T>::operator[](index_type i)				{ __ASSERT(i < _size); return data[i]; }
template<typename T> inline T &Grid3x4<T>::operator()(index_type i)				{ __ASSERT(i < _size); return data[i]; }

template<typename T> inline const T Grid3x4<T>::operator()(index_type i, index_type j, index_type k, index_type l) const { 
	__ASSERT(i < nX && j < nY && k < nZ && l < 4);
	//return data[4*((i*nY+j)*nZ+k)+l]; 
	return data[l+4*(k+nZ*(j+i*nY))]; 
}

template<typename T> inline T &Grid3x4<T>::operator()(index_type i, index_type j, index_type k, index_type l) { 
	__ASSERT(i < nX && j < nY && k < nZ && l < 4);
	//return data[4*((i*nY+j)*nZ+k)+l]; 
	return data[l+4*(k+nZ*(j+i*nY))]; 
}

template<typename T> inline const T Grid3x4<T>::get(index_type i, index_type j, index_type k, index_type l) const { 
	__ASSERT(i < nX && j < nY && k < nZ && l < 4);
	//return data[4*((i*nY+j)*nZ+k)+l]; 
	return data[l+4*(k+nZ*(j+i*nY))]; 
}

//////////////////////////////////////////////////////////////////////////
// Assignment and set
template<typename T> void Grid3x4<T>::zero()					{ VectorOperation<Grid3x4<T>>().zero(*this);		}

template<typename T> void Grid3x4<T>::set(const Grid3x4 &grid)	{ VectorOperation<Grid3x4<T>>().set(*this, grid);	}
template<typename T> void Grid3x4<T>::set(const T &val)			{ VectorOperation<Grid3x4<T>>().set(*this, val);	}

template<typename T> const Grid3x4<T> &Grid3x4<T>::operator=(const Grid3x4 &grid)	{ set(grid); return (*this); }

//////////////////////////////////////////////////////////////////////////
// Manipulation	
template <typename T> const T Grid3x4<T>::infNorm() const {
	VectorFunctor_InfNorm<Grid3x4> func(*this);
	func(blocked_range<typename Grid3x4<T>::index_type>(0,_size));
	return func.getValue();	
}

template <typename T> const T Grid3x4<T>::l2Norm() const {
	VectorFunctor_L2Norm<Grid3x4> func(*this);
	func(blocked_range<typename Grid3x4<T>::index_type>(0,_size));
	return func.getValue();	
}

template<typename T> const Grid3x4<T> &Grid3x4<T>::operator*=(const T c) { 
	VectorFunctor_Scale<Grid3x4,T>(*this, c)(blocked_range<index_type>(0,_size));
	return (*this); 
}	

#undef __ASSERT

};	// namespace arithmetic

#endif