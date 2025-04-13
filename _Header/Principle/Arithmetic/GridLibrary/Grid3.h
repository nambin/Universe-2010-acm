#ifndef _ARITHMETIC_GRID3_H_
#define _ARITHMETIC_GRID3_H_

#include "../../../Basis/Basis/Environments.h"
#include "../../../Basis/Basis/BasicLibrary.h"
#include "../../../Basis/Basis/Functors.h"
#include "../../Geometry/GeometricPrimitive/AABB.h"
#include "../Vector/Vector3.h"
#include "GridLibrary.h"
#include "Grid3_Functors.h"
#include "ArrayFunctor.h"
#include "GridBoundaryCondition.h"

namespace arithmetic {

using namespace basis;
using arithmetic::Vector3;

#ifdef _DEBUG
#define __ASSERT(c) assert(c)
#else
#define __ASSERT(c) {}
#endif

template <typename T, typename R> class Grid3 {
public : 
	typedef grid_index_type		index_type;
	typedef index_type*			index_pointer ;
	typedef const index_type*	index_const_pointer ;
	typedef index_type&			index_reference ;
	typedef const index_type&	index_const_reference ;

	typedef T					value_type;
	typedef value_type*			value_pointer ;       
	typedef const value_type*	value_const_pointer ; 
	typedef value_type&			value_reference ;     
	typedef const value_type&	value_const_reference;

	typedef R					real_type;

	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Grid3();
	Grid3(const Grid3 &grid3);
	Grid3(index_type nX, index_type nY, index_type nZ, R dx, R dy, R dz, Vector3<R> offset = Vector3<R>(0.0));
	~Grid3();

	void init(index_type nX, index_type nY, index_type nZ, R dx, R dy, R dz, const Vector3<R> &offset = Vector3<R>(0.0));
	template<typename S, typename R2> void resize(const Grid3<S,R2>&);
	void release();

	void set_offset(const Vector3<R>&);
	void set_origin(const Vector3<R>&);

	index_type size() const { return _size; }

	template<int DIM> const bool is_valid_range(index_type) const;
	template<int DIM> const bool is_valid_range_for_interp(index_type) const;

	const bool is_empty() const;

	//////////////////////////////////////////////////////////////////////////
	// Index operators
	inline const index_type idx(index_type i, index_type j, index_type k) const;
	inline const index_type idx(const Vector3<index_type>& idxvec) const;

	inline const T	operator[](index_type i) const ;
	inline T&		operator[](index_type i)		;
	inline const T	operator()(index_type i) const ;
	inline T&		operator()(index_type i)		;
	inline const T	operator()(index_type i, index_type j, index_type k) const	;
	inline T&		operator()(index_type i, index_type j, index_type k)		;
	inline const T	operator()(const Vector3<index_type>&) const;
	inline T&		operator()(const Vector3<index_type>&)		;

	inline const T	get(index_type i, index_type j, index_type k) const;
	inline const T	get(const Vector3<index_type>&) const;

	//////////////////////////////////////////////////////////////////////////
	// Getter	
	inline const T get(const Vector3<R> &_pos) const;
	inline const T operator()(const Vector3<R> &_pos) const;

	inline const T	get_nearest(const Vector3<R>&) const;

	inline const Vector3<R> pos(index_type i, index_type j, index_type k) const;
	inline const Vector3<R> pos(const index_type idxvec[3]) const;
	inline const Vector3<R> pos(const Vector3<index_type>&) const;

	inline void getMaxAndMin(T &_max, T &_min, const Vector3<R> &_pos) const;

	//////////////////////////////////////////////////////////////////////////
	// Get Pointer
	T* getptr() const { return data; }
	
	//////////////////////////////////////////////////////////////////////////
	// Assignment and set
	void zero();

	template<typename T1, typename R1> void set(const Grid3<T1,R1> &grid);
	void set(const Grid3 &grid);
	void set(const R* da);
	void set(const T &val);

	void set_boundary(const T& val);

	const Grid3& operator=(const Grid3 &grid);

	//////////////////////////////////////////////////////////////////////////
	// Manipulation	
	const Grid3& setNeumannBoundary();
	const T infNorm() const;
	const T l2Norm() const;
	const T sum() const;

	//////////////////////////////////////////////////////////////////////////
	// In place arithmetic
	const Grid3& operator+=(const Grid3 &grid);
	const Grid3& operator-=(const Grid3 &grid);
	const Grid3& operator*=(const R c);
	
	void scale(R c) { (*this) *= c;						}
	void negate()	{ (*this) *= static_cast<R>(-1.0);	}

	//////////////////////////////////////////////////////////////////////////
	// Utility
	inline const geometryspace::AABB<R> getAABB() const;

	//////////////////////////////////////////////////////////////////////////
	// Friend Methods
	template <typename T, typename R> friend void swap(Grid3<T,R> &a, Grid3<T,R> &b);	
	template <typename T, typename R> friend const T dot(const Grid3<T,R> &a, const Grid3<T,R> &b);	

public : 
	void get_idx_normalized_coord(index_type &i, index_type &j, index_type& k, R& fx, R& fy, R& fz, const Vector3<R>& _pos) const;

	inline void getIdx(index_type &i, index_type &j, index_type &k, const Vector3<R>& _pos) const;
	inline void getClosestIdx(index_type &i, index_type &j, index_type& k, const Vector3<R>& _pos) const;
	inline const T trilerp(index_type i, index_type j, index_type k, R fx, R fy, R fz) const;

	inline void getIdx(Vector3<index_type>& idx, const Vector3<R>& _pos) const				{ getIdx(idx[0], idx[1], idx[2], _pos);					}
	inline void getClosestIdx(Vector3<index_type>& idx, const Vector3<R>& _pos) const		{ getClosestIdx(idx[0], idx[1], idx[2], _pos);			}
	inline const T trilerp(const Vector3<index_type>& idx, R fx, R fy, R fz) const			{ return trilerp(idx[0], idx[1], idx[2], fx, fy, fz);	}

	const bool exportBinary(std::ofstream&) const;
	const bool importBinary(std::ifstream&);

public : 
	index_type	dim[3];
	index_type	nX, nY, nZ;
	index_type	_size;

	T *			data;

	R dx, dy, dz;
	R dxInv, dyInv, dzInv;
	R dxInvHalf, dyInvHalf, dzInvHalf;

	Vector3<R> offset, origin;
};

template<typename T, typename R> 
const bool Grid3<T,R>::exportBinary(std::ofstream& ofs) const {
	ofs.write((char*)(dim), sizeof(index_type) * 3);
	ofs.write((char*)(&nX), sizeof(index_type));
	ofs.write((char*)(&nY), sizeof(index_type));
	ofs.write((char*)(&nZ), sizeof(index_type));
	ofs.write((char*)(&_size), sizeof(index_type));

	ofs.write((char*)(data), sizeof(T) * _size);

	ofs.write((char*)(&dx), sizeof(R));			ofs.write((char*)(&dy), sizeof(R));			ofs.write((char*)(&dz), sizeof(R));
	ofs.write((char*)(&dxInv), sizeof(R));		ofs.write((char*)(&dyInv), sizeof(R));		ofs.write((char*)(&dzInv), sizeof(R));
	ofs.write((char*)(&dxInvHalf), sizeof(R));	ofs.write((char*)(&dyInvHalf), sizeof(R));	ofs.write((char*)(&dzInvHalf), sizeof(R));

	ofs.write((char*)(&offset), sizeof(R)*3);
	ofs.write((char*)(&origin), sizeof(R)*3);

	return true;
}

template<typename T, typename R> 
const bool Grid3<T,R>::importBinary(std::ifstream& ifs) {
	ifs.read((char*)(dim), sizeof(index_type) * 3);
	ifs.read((char*)(&nX), sizeof(index_type));
	ifs.read((char*)(&nY), sizeof(index_type));
	ifs.read((char*)(&nZ), sizeof(index_type));
	ifs.read((char*)(&_size), sizeof(index_type));

	if(data) delete data; data = new T[_size];
	ifs.read((char*)(data), sizeof(T) * _size);

	ifs.read((char*)(&dx), sizeof(R));			ifs.read((char*)(&dy), sizeof(R));			ifs.read((char*)(&dz), sizeof(R));
	ifs.read((char*)(&dxInv), sizeof(R));		ifs.read((char*)(&dyInv), sizeof(R));		ifs.read((char*)(&dzInv), sizeof(R));
	ifs.read((char*)(&dxInvHalf), sizeof(R));	ifs.read((char*)(&dyInvHalf), sizeof(R));	ifs.read((char*)(&dzInvHalf), sizeof(R));

	ifs.read((char*)(&offset), sizeof(R)*3);
	ifs.read((char*)(&origin), sizeof(R)*3);

	return true;
}

//////////////////////////////////////////////////////////////////////////
// Constructor/Destructor
template<typename T, typename R> 
Grid3<T,R>::Grid3() 
: nX(0), nY(0), nZ(0), _size(0), data(NULL), dx(0.0), dy(0.0), dz(0.0), dxInv(0.0), dyInv(0.0), dzInv(0.0)
, dxInvHalf(0.0), dyInvHalf(0.0), dzInvHalf(0.0), offset(0), origin(0)
{ dim[0] = dim[1] = dim[2] = 0; }

template<typename T, typename R> 
Grid3<T,R>::Grid3(const Grid3<T,R> &grid3)
: nX(0), nY(0), nZ(0), _size(0), data(NULL), dx(0.0), dy(0.0), dz(0.0), dxInv(0.0), dyInv(0.0), dzInv(0.0)
, dxInvHalf(0.0), dyInvHalf(0.0), dzInvHalf(0.0), offset(0), origin(0) 
{
	dim[0] = dim[1] = dim[2] = 0;
	init(grid3.nX, grid3.nY, grid3.nZ, grid3.dx, grid3.dy, grid3.dz, grid3.offset);
	set(grid3);
}

template<typename T, typename R> 
Grid3<T,R>::Grid3(index_type _nX, index_type _nY, index_type _nZ, R _dx, R _dy, R _dz, Vector3<R> _offset) 
: nX(0), nY(0), nZ(0), _size(0), data(NULL), dx(0.0), dy(0.0), dz(0.0), dxInv(0.0), dyInv(0.0), dzInv(0.0)
, dxInvHalf(0.0), dyInvHalf(0.0), dzInvHalf(0.0), offset(0), origin(0)
{ dim[0] = dim[1] = dim[2] = 0; init(_nX, _nY, _nZ, _dx, _dy, _dz, _offset); }

template<typename T, typename R> 
Grid3<T,R>::~Grid3()	
{ release(); }

template<typename T, typename R> void Grid3<T,R>::init(index_type _nX, index_type _nY, index_type _nZ, R _dx, R _dy, R _dz, const Vector3<R> &_offset)	{
	release();

	dim[0] = _nX; dim[1] = _nY; dim[2] = _nZ;
	nX = _nX; nY = _nY; nZ = _nZ;
	_size = nX * nY * nZ;

	data = new T[_size];
	//data = (T*)memoryspace::AllocAligned(_size*sizeof(T));

	dx = _dx; dy = _dy;	dz = _dz;

	if(IsAlmostZero(dx) || IsAlmostZero(dy) || IsAlmostZero(dz)) { dxInv = 0; dyInv = 0; dzInv = 0; }
	else { dxInv = R(1.0)/dx;	dyInv = R(1.0)/dy; dzInv = R(1.0)/dz; dxInvHalf = R(0.5) * dxInv; dyInvHalf = R(0.5) * dyInv; dzInvHalf = R(0.5) * dzInv; }

	offset = _offset;
	origin = Vector3<R>(offset[0] * dx, offset[1] * dy, offset[2] * dz);
}

template<typename T, typename R> template<typename S, typename R2> void Grid3<T,R>::resize(const Grid3<S,R2>& grid) {
	init(grid.nX, grid.nY, grid.nZ, grid.dx, grid.dy, grid.dz, Vector3<R>(grid.offset));
}

template<typename T, typename R> void Grid3<T,R>::release()	{
	if(data) delete[] data; data = NULL;
	//if(data) memoryspace::FreeAligned(data); data = NULL;
	nX = nY = nZ = _size = 0;
}

template<typename T, typename R> void Grid3<T,R>::set_offset(const Vector3<R>& _offset) {
	offset = _offset;
	origin = Vector3<R>(offset[0] * dx, offset[1] * dy, offset[2] * dz);
}

template<typename T, typename R> void Grid3<T,R>::set_origin(const Vector3<R>& _origin) {
	if(IsAlmostZero(dx) || IsAlmostZero(dy) || IsAlmostZero(dz))
		basis::throwError("Should Initialize Grid3 first\n");

	origin = _origin;
	offset = Vector3<R>(origin[0] / dx, origin[1] / dy, origin[2] / dz);
}

template<typename T, typename R> template<int DIM> const bool Grid3<T,R>::is_valid_range(index_type i) const {
	if(0 <= i && i < dim[DIM])	return true;
	else						return false;
}

template<typename T, typename R> template<int DIM> const bool Grid3<T,R>::is_valid_range_for_interp(index_type i) const {
	if(0 <= i && i < dim[DIM]-1)	return true;
	else							return false;
}

template<typename T, typename R> const bool Grid3<T,R>::is_empty() const {
	if(nX == 0 || nY == 0 || nZ == 0)
		return true;
	return false;
}

//////////////////////////////////////////////////////////////////////////
// Index operators
template<typename T, typename R> inline const T Grid3<T,R>::operator[](index_type i) const	{ __ASSERT(i < _size); return data[i];	}
template<typename T, typename R> inline const T Grid3<T,R>::operator()(index_type i) const	{ __ASSERT(i < _size); return data[i];	}
template<typename T, typename R> inline T& Grid3<T,R>::operator[](index_type i)				{ __ASSERT(i < _size); return data[i];	}
template<typename T, typename R> inline T& Grid3<T,R>::operator()(index_type i)				{ __ASSERT(i < _size); return data[i];	}

template<typename T, typename R> inline const typename Grid3<T,R>::index_type Grid3<T,R>::idx(index_type i, index_type j, index_type k) const { 
	__ASSERT(i < nX && j < nY && k < nZ);
	return ((i*nY+j)*nZ+k);							
}

template<typename T, typename R> inline const typename Grid3<T,R>::index_type Grid3<T,R>::idx(const Vector3<index_type>& idx) const { 
	__ASSERT(idx[0] < nX && idx[1] < nY && idx[2] < nZ);
	return ((idx[0]*nY+idx[1])*nZ+idx[2]);	
}

template<typename T, typename R> inline const T Grid3<T,R>::operator()(index_type i, index_type j, index_type k) const 	{ 
	__ASSERT(i < nX && j < nY && k < nZ);
	return data[(i*nY+j)*nZ+k]; 
}

template<typename T, typename R> inline T& Grid3<T,R>::operator()(index_type i, index_type j, index_type k)	{ 
	__ASSERT(i < nX && j < nY && k < nZ);
	return data[(i*nY+j)*nZ+k]; 
}

template<typename T, typename R> inline const T Grid3<T,R>::operator()(const Vector3<index_type>& idx) const 	{ 
	__ASSERT(idx[0] < nX && idx[1] < nY && idx[2] < nZ);
	return data[(idx[0]*nY+idx[1])*nZ+idx[2]]; 
}

template<typename T, typename R> inline T& Grid3<T,R>::operator()(const Vector3<index_type>& idx)	{ 
	__ASSERT(idx[0] < nX && idx[1] < nY && idx[2] < nZ);
	return data[(idx[0]*nY+idx[1])*nZ+idx[2]]; 
}

template<typename T, typename R> inline const T Grid3<T,R>::get(index_type i, index_type j, index_type k) const	{ 
	__ASSERT(i < nX && j < nY && k < nZ);
	return data[(i*nY+j)*nZ+k]; 
}

template<typename T, typename R> inline const T Grid3<T,R>::get(const Vector3<index_type>& idx)	const { 
	__ASSERT(idx[0] < nX && idx[1] < nY && idx[2] < nZ);
	return data[(idx[0]*nY+idx[1])*nZ+idx[2]]; 
}

//////////////////////////////////////////////////////////////////////////
// Getter
template<typename T, typename R> inline const T Grid3<T,R>::get(const Vector3<R> &_pos) const {
	return (*this)(_pos);
}

template<typename T, typename R> inline const T Grid3<T,R>::operator()(const Vector3<R> &_pos) const { 
	index_type i,j,k;
	getIdx(i,j,k,_pos);
	const Vector3<R> corner = pos(i,j,k);
	const R fx = (_pos[0] - corner[0]) * dxInv;
	const R fy = (_pos[1] - corner[1]) * dyInv;
	const R fz = (_pos[2] - corner[2]) * dzInv;
	return trilerp(i,j,k,fx,fy,fz);
}

template<typename T, typename R> inline const T Grid3<T,R>::get_nearest(const Vector3<R>& _pos) const {
	index_type i,j,k;
	getIdx(i,j,k,_pos);
	const Vector3<R> corner = pos(i,j,k);
	const R fx = (_pos[0] - corner[0]) * dxInv;
	const R fy = (_pos[1] - corner[1]) * dyInv;
	const R fz = (_pos[2] - corner[2]) * dzInv;

	index_type i_icr = (fx < 0.5) ? 0 : 1;
	index_type j_icr = (fy < 0.5) ? 0 : 1;
	index_type k_icr = (fz < 0.5) ? 0 : 1;

	return (*this)(i+i_icr,j+j_icr,k+k_icr);
}

template<typename T, typename R> inline const Vector3<R> Grid3<T,R>::pos(index_type i, index_type j, index_type k) const {
	return Vector3<R>(i*dx+origin[0], j*dy+origin[1], k*dz+origin[2]); 
}

template<typename T, typename R> inline const Vector3<R> Grid3<T,R>::pos(const index_type i[3]) const {
	return Vector3<R>(i[0]*dx+origin[0], i[1]*dy+origin[1], i[2]*dz+origin[2]); 
}

template<typename T, typename R> inline const Vector3<R> Grid3<T,R>::pos(const Vector3<index_type>& idx) const {
	return Vector3<R>(idx[0]*dx+origin[0], idx[1]*dy+origin[1], idx[2]*dz+origin[2]); 
}

template<typename T, typename R> inline void Grid3<T,R>::getMaxAndMin(T &_max, T &_min, const Vector3<R> &_pos) const {
	index_type i,j,k;
	getIdx(i,j,k,_pos);
	_max = MAX(	MAX(	MAX((*this)(i,j,k),(*this)(i+1,j,k)), 
						MAX((*this)(i,j+1,k),(*this)(i+1,j+1,k)) ), 
				MAX(	MAX((*this)(i,j,k+1),(*this)(i+1,j,k+1)), 
						MAX((*this)(i,j+1,k+1),(*this)(i+1,j+1,k+1)) ) );
	_min = MIN(	MIN(	MIN((*this)(i,j,k),(*this)(i+1,j,k)), 
						MIN((*this)(i,j+1,k),(*this)(i+1,j+1,k)) ),
				MIN(	MIN((*this)(i,j,k+1),(*this)(i+1,j,k+1)), 
						MIN((*this)(i,j+1,k+1),(*this)(i+1,j+1,k+1)) ) );
}

//////////////////////////////////////////////////////////////////////////
// Assignment and set
template<typename T, typename R> void Grid3<T,R>::zero()						{ VectorOperation<Grid3<T,R>>().zero(*this); }

template<typename T, typename R> template<typename T1, typename R1> void Grid3<T,R>::set(const Grid3<T1,R1> &grid) {
	if(nX != grid.nX || nY != grid.nY || nZ != grid.nZ)
		throwError("Error in Grid3<T,R>::void set(const Grid3<T1,R1> &grid)");

	for(index_type i=0;i<nX;++i)
		for(index_type j=0;j<nY;++j)
			for(index_type k=0;k<nZ;++k)
				(*this)(i,j,k) = grid(i,j,k);
}

template<typename T, typename R> void Grid3<T,R>::set(const Grid3<T,R> &grid)	{ VectorOperation<Grid3<T,R>>().set(*this, grid);	}
template<typename T, typename R> void Grid3<T,R>::set(const R* da)				{ VectorOperation<Grid3<T,R>>().set(*this, da);		}
template<typename T, typename R> void Grid3<T,R>::set(const T &val)				{ VectorOperation<Grid3<T,R>>().set(*this, val);	}

template<typename T, typename R> void Grid3<T,R>::set_boundary(const T& val) {

#ifdef USE_OPENMP
	int my_end = int(nX);
	int i;

#pragma omp parallel for
#else
	index_type my_end = nX;
	index_type i;
#endif
	for( i = 0; i < my_end; ++i ) 	
		for(index_type j=0;j<nY;j++)
			(*this)(i,j,0) = (*this)(i,j,nZ-1) = true;

#ifdef USE_OPENMP
	my_end = int(nX);	
#pragma omp parallel for
#else
	my_end = nX;	
#endif
	for( i = 0; i < my_end; ++i ) 
		for(index_type k=0;k<nZ;k++)
			(*this)(i,0,k) = (*this)(i,nY-1,k) = true;

#ifdef USE_OPENMP
	my_end = int(nY);	
	int j;
#pragma omp parallel for
#else
	my_end = nY;	
	index_type j;
#endif
	for( j = 0; j < my_end; ++j ) 
		for(index_type k=0;k<nZ;k++)
			(*this)(0,j,k) = (*this)(nX-1,j,k) = true;
}

template<typename T, typename R> const Grid3<T,R> &Grid3<T,R>::operator=(const Grid3<T,R> &grid) { set(grid); return (*this); }

//////////////////////////////////////////////////////////////////////////
// Manipulation	
template <typename T, typename R> const Grid3<T,R>& Grid3<T,R>::setNeumannBoundary()	{
	GridBoundaryCondition_Neumann<T,R> neumann; neumann.enforce(*this);
	return (*this);
}

template <typename T, typename R> const T Grid3<T,R>::infNorm()	const {
	VectorFunctor_InfNorm<Grid3> func(*this);
	func(blocked_range<index_type>(0,_size));
	return func.getValue();	
}

template <typename T, typename R> const T Grid3<T,R>::l2Norm() const {
	VectorFunctor_L2Norm<Grid3> func(*this);
	func(blocked_range<index_type>(0,_size));
	return func.getValue();	
}

template <typename T, typename R> const T Grid3<T,R>::sum() const {
	VectorFunctor_Sum<Grid3> func(*this);
	func(blocked_range<index_type>(0,_size));
	return func.getValue();	
}

//////////////////////////////////////////////////////////////////////////
// In place arithmetic
template<typename T, typename R> const Grid3<T,R> &Grid3<T,R>::operator+=(const Grid3 &grid) { VectorFunctor_Icr<Grid3>(*this, g)(blocked_range<index_type>(0,_size)); return (*this);	}
template<typename T, typename R> const Grid3<T,R> &Grid3<T,R>::operator-=(const Grid3 &grid) { VectorFunctor_Dcr<Grid3>(*this, g)(blocked_range<index_type>(0,_size)); return (*this);	}
template<typename T, typename R> const Grid3<T,R> &Grid3<T,R>::operator*=(const R c)		 { VectorFunctor_Scale<Grid3,R>(*this, c)(blocked_range<index_type>(0,_size));return (*this); }	

template<typename T, typename R> inline const geometryspace::AABB<R> Grid3<T,R>::getAABB() const {
	return geometryspace::AABB<R>(origin, origin + Vector3<R>(dx * (nX-1), dy * (nY-1), dz * (nZ-1)));
}

//////////////////////////////////////////////////////////////////////////
// Private Functions
template<typename T, typename R> inline void Grid3<T,R>::getIdx(index_type &i, index_type &j, index_type &k, const Vector3<R>& _pos) const {
	i = (index_type)((_pos[0] - origin[0]) * dxInv);
	j = (index_type)((_pos[1] - origin[1]) * dyInv);
	k = (index_type)((_pos[2] - origin[2]) * dzInv);

	if(i == nX-1) --i;	if(j == nY-1) --j; if(k == nZ-1) --k;
}

template<typename T, typename R> void Grid3<T,R>::get_idx_normalized_coord(index_type &i, index_type &j, index_type& k, R& fx, R& fy, R& fz, const Vector3<R>& _pos) const {
	getIdx(i,j,k,_pos);
	Vector3<R> _corner = pos(i,j,k);
	fx = (_pos[0] - _corner[0]) * dxInv;
	fy = (_pos[1] - _corner[1]) * dyInv;
	fz = (_pos[2] - _corner[2]) * dzInv;
}

template<typename T, typename R> inline void Grid3<T,R>::getClosestIdx(index_type &i, index_type &j, index_type& k, const Vector3<R>& _pos) const {
	R ir = ((_pos[0] - origin[0]) * dxInv);
	R jr = ((_pos[1] - origin[1]) * dyInv);
	R kr = ((_pos[2] - origin[2]) * dzInv);

	i = (index_type)ir; j = (index_type)jr; k = (index_type)kr;

	if(ir - i > 0.5f) ++i;
	if(jr - j > 0.5f) ++j;
	if(kr - k > 0.5f) ++k;

	if(i == nX-1) --i;	else if(i == 0) ++i;
	if(j == nY-1) --j;	else if(j == 0)	++j;
	if(k == nZ-1) --k;	else if(k == 0)	++k;
}

template<typename T, typename R> inline const T Grid3<T,R>::trilerp(index_type i, index_type j, index_type k, R fx, R fy, R fz) const {
	return (1-fx) * ((1-fy) * ((1-fz) * (*this)(i,j,k)     + fz * (*this)(i,j,k+1))   
					  + fy  * ((1-fz) * (*this)(i,j+1,k)   + fz * (*this)(i,j+1,k+1)))
			+ fx  * ((1-fy) * ((1-fz) * (*this)(i+1,j,k)   + fz * (*this)(i+1,j,k+1)) 
					  + fy  * ((1-fz) * (*this)(i+1,j+1,k) + fz * (*this)(i+1,j+1,k+1)));
}

template <typename T, typename R> inline void swap(Grid3<T,R> &a, Grid3<T,R> &b)	{
	T *tmp = a.data;	a.data = b.data;	b.data = tmp;	
}

template <typename T, typename R> const T dot(const Grid3<T,R> &a, const Grid3<T,R> &b)	{
	__ASSERT(a._size == b._size);
	ArrayDotFunctor<T> dotfunc(a.data, b.data);
	dotfunc(basis::blocked_range<typename Grid3<T,R>::index_type>(0, a._size));
	return dotfunc.getValue();
}

#undef __ASSERT

};		// namespace arithmetic

#endif