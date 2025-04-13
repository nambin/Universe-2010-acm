#ifndef _ARITHMETIC_GRID2x3_H_
#define _ARITHMETIC_GRID2x3_H_

#include <fstream>
#include "../../../Basis/Basis/Environments.h"
#include "../../../Basis/Basis/BasicLibrary.h"
#include "Grid2.h"

namespace arithmetic {

#ifdef _DEBUG
#define __ASSERT(c) assert(c)
#else
#define __ASSERT(c) {}
#endif

template <typename T> class Grid2x3	{
public : 
	typedef typename Grid2<T,T>::index_type	index_type;
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
	Grid2x3();
	Grid2x3(index_type _nX, index_type _nY);
	~Grid2x3();

	void init(index_type _nX, index_type _nY);
	void set(index_type _nX, index_type _nY);
	void release();

	index_type size() const { return _size; }

	//////////////////////////////////////////////////////////////////////////
	// Index operators
	inline const T	operator[](index_type i) const;
	inline T&		operator[](index_type i)		;
	inline const T	operator()(index_type i) const;
	inline T&		operator()(index_type i)		;
	inline const T	operator()(index_type i, index_type j, index_type k) const;
	inline T&		operator()(index_type i, index_type j, index_type k)		;

	inline const T	get(index_type i, index_type j, index_type k) const;

	//////////////////////////////////////////////////////////////////////////
	// Get Pointer
	T* getptr() const { return data; }
	
	//////////////////////////////////////////////////////////////////////////
	// Assignment and set
	void zero();

	void set(const Grid2x3 &grid);
	void set(const T &val);

	const Grid2x3 &operator=(const Grid2x3 &grid);

	//////////////////////////////////////////////////////////////////////////
	// Manipulation		
	const T l2Norm() const;
	const T infNorm() const;

	const Grid2x3 &operator*=(const T c);

	const bool exportBinary(std::ofstream&) const;
	const bool importBinary(std::ifstream&);

public : 
	index_type	nX, nY;
	index_type	_size;
	T *			data;	
};

template<typename T>
const bool Grid2x3<T>::exportBinary(std::ofstream& ofs) const {
	ofs.write((char*)(&nX), sizeof(index_type));
	ofs.write((char*)(&nY), sizeof(index_type));
	ofs.write((char*)(&_size), sizeof(index_type));
	ofs.write((char*)(data), sizeof(T) * _size);

	return true;
}

template<typename T>
const bool Grid2x3<T>::importBinary(std::ifstream& ifs) {
	ifs.read((char*)(&nX), sizeof(index_type));
	ifs.read((char*)(&nY), sizeof(index_type));
	ifs.read((char*)(&_size), sizeof(index_type));

	if(data) delete data; data = new T[_size];
	ifs.read((char*)(data), sizeof(T) * _size);

	return true;
}

//////////////////////////////////////////////////////////////////////////
// Constructor/Destructor
template<typename T> Grid2x3<T>::Grid2x3()
: nX(0), nY(0), _size(0), data(NULL)
{}

template<typename T> Grid2x3<T>::Grid2x3(index_type _nX, index_type _nY)
: nX(0), nY(0), _size(0), data(NULL) { 
	set(_nX, _nY); 
}

template<typename T> Grid2x3<T>::~Grid2x3() { 
	release(); 
}

template<typename T> void Grid2x3<T>::init(index_type _nX, index_type _nY) {
	set(_nX, _nY);
}

template<typename T> void Grid2x3<T>::set(index_type _nX, index_type _nY)	{
	nX = _nX; nY = _nY;

	release();
	nX = _nX;
	nY = _nY;
	_size = nX * nY * 3;
	data = new T[_size];
	zero();
}

template<typename T> void Grid2x3<T>::release()	{
	if(data) delete[] data; data = NULL;
	nX = nY = _size = 0;
}

//////////////////////////////////////////////////////////////////////////
// Index operators
template<typename T> inline T &Grid2x3<T>::operator[](index_type i)				{ __ASSERT(i < _size); return data[i];	}
template<typename T> inline T &Grid2x3<T>::operator()(index_type i)				{ __ASSERT(i < _size); return data[i];	}
template<typename T> inline const T Grid2x3<T>::operator[](index_type i) const	{ __ASSERT(i < _size); return data[i];	}
template<typename T> inline const T Grid2x3<T>::operator()(index_type i) const	{ __ASSERT(i < _size); return data[i];	}

template<typename T> inline const T Grid2x3<T>::operator()(index_type i, index_type j, index_type k) const { 
	__ASSERT(i < nX && j < nY && k < 3);
	return data[3*(i*nY+j)+k];  
}

template<typename T> inline T &Grid2x3<T>::operator()(index_type i, index_type j, index_type k) { 
	__ASSERT(i < nX && j < nY && k < 3);
	return data[3*(i*nY+j)+k];
}

template<typename T> inline const T Grid2x3<T>::get(index_type i, index_type j, index_type k) const { 
	assert(i < nX && j < nY && k < 3);
	return data[3*(i*nY+j)+k];
}

//////////////////////////////////////////////////////////////////////////
// Assignment and set
template<typename T> void Grid2x3<T>::zero()						{ VectorOperation<Grid2x3<T>>().zero(*this);		}

template<typename T> void Grid2x3<T>::set(const Grid2x3 &grid)		{ VectorOperation<Grid2x3<T>>().set(*this, grid);	}
template<typename T> void Grid2x3<T>::set(const T &val)				{ VectorOperation<Grid2x3<T>>().set(*this, val);	}

template<typename T> const Grid2x3<T> &Grid2x3<T>::operator=(const Grid2x3 &grid)	{ set(grid); return grid; }

//////////////////////////////////////////////////////////////////////////
// Manipulation	
template <typename T> const T Grid2x3<T>::l2Norm() const {
	VectorFunctor_L2Norm<Grid2x3> func(*this);
	func(blocked_range<typename Grid2x3<T>::index_type>(0,_size));
	return func.getValue();	
}

template <typename T> const T Grid2x3<T>::infNorm() const {
	VectorFunctor_InfNorm<Grid2x3> func(*this);
	func(blocked_range<typename Grid2x3<T>::index_type>(0,_size));
	return func.getValue();	
}

template<typename T> const Grid2x3<T> &Grid2x3<T>::operator*=(const T c) { 
	VectorFunctor_Scale<Grid2x3,T>(*this, c)(blocked_range<index_type>(0,_size));
	return (*this); 
}	

#undef __ASSERT

};	// namespace arithmetic

#endif