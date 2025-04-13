#ifndef _ARITHMETIC_GRID2_H_
#define _ARITHMETIC_GRID2_H_

#include "../../../Basis/Basis/Environments.h"
#include "../../../Basis/Basis/BasicLibrary.h"
#include "../../../Basis/Basis/Functors.h"
#include "../../../Basis/Color/RGBA.h"
#include "../../Geometry/GeometricPrimitive/AABB2D.h"
#include "../Vector/Vector2.h"
#include "../Vector/VectorOperation.h"
#include "GridLibrary.h"
#include "Grid2_Functors.h"
#include "ArrayFunctor.h"
#include "GridBoundaryCondition.h"

namespace arithmetic {

using namespace basis;
using arithmetic::Vector2;
using basis::blocked_range;
using basis::blocked_range2d;

#ifdef _DEBUG
#define __ASSERT(c) assert(c)
#else
#define __ASSERT(c) {}
#endif

template <typename T, typename R> class Grid2 {
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
	Grid2();
	Grid2(const Grid2 &grid2);
	Grid2(index_type nX, index_type nY, R dx, R dy, Vector2<R> offset = Vector2<R>(0.0));
	~Grid2();

	void init(index_type nX, index_type nY, R dx, R dy, const Vector2<R> &offset = Vector2<R>(0.0));
	template<typename S, typename R2> void resize(const Grid2<S,R2>&);
	void release();	

	void set_offset(const Vector2<R>&);
	void set_origin(const Vector2<R>&);

	index_type size() const { return _size; }

	template<int DIM> const bool is_valid_range(index_type) const;
	template<int DIM> const bool is_valid_range_for_interp(index_type) const;

	const bool is_empty() const;

	//////////////////////////////////////////////////////////////////////////
	// Index operators
	inline const index_type idx(index_type i, index_type j) const			{ return (i*nY+j); }
	inline const index_type idx(const Vector2<index_type>& idxvec) const	{ return (idxvec[0]*nY+idxvec[1]); }

	inline const T	operator[](index_type i) const	;
	inline T&		operator[](index_type i)		;
	inline const T	operator()(index_type i) const	;
	inline T&		operator()(index_type i)		;
	inline const T	operator()(index_type i, index_type j) const;
	inline T&		operator()(index_type i, index_type j)		;
	inline const T	operator()(const Vector2<index_type>& idxvec) const	;
	inline T&		operator()(const Vector2<index_type>& idxvec)		;

	inline const T	get(index_type i, index_type j) const;
	inline const T	get(const Vector2<index_type>& idxvec) const;
	
	//////////////////////////////////////////////////////////////////////////
	// Getter	
	inline const T	get(const Vector2<R> &_pos) const;
	inline const T	operator()(const Vector2<R> &_pos) const;

	inline const T	get_nearest(const Vector2<R>&) const;

	inline const Vector2<R> pos(index_type i, index_type j) const;
	inline const Vector2<R> pos(const index_type idxvec[2]) const;
	inline const Vector2<R> pos(const Vector2<index_type>& idxvec) const;

	inline void getMaxAndMin(T &_max, T &_min, const Vector2<R> &_pos) const;	

	//////////////////////////////////////////////////////////////////////////
	// Get Pointer
	T* getptr() const { return data; }

	//////////////////////////////////////////////////////////////////////////
	// Assignment and set
	void zero();

	template<typename T1, typename R1> void set(const Grid2<T1,R1> &grid);
	void set(const Grid2 &grid);
	void set(const T val);

	void set_boundary(const T val);

	const Grid2 &operator=(const Grid2 &grid);

	//////////////////////////////////////////////////////////////////////////
	// Manipulation	
	void setNeumannBoundary();
	const T infNorm() const;
	const T l2Norm() const;
	const T sum() const;

	//////////////////////////////////////////////////////////////////////////
	// In place arithmetic
	const Grid2 &operator+=(const Grid2 &grid);
	const Grid2 &operator-=(const Grid2 &grid);
	const Grid2 &operator*=(const R c);

	void scale(R c) { (*this) *= c;						}
	void negate()	{ (*this) *= static_cast<R>(-1);	}

	//////////////////////////////////////////////////////////////////////////
	// Utility
	inline const geometryspace::AABB2D<R> getAABB() const;

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	void glDrawGrid(const RGBA &color, GLfloat width) const;

	//////////////////////////////////////////////////////////////////////////
	// IO
	void save_as_matlab(const char* const filename, const char* const grid_name = "phi") const;

	//////////////////////////////////////////////////////////////////////////
	// Friend Methods
	template <typename T, typename R> friend void swap(Grid2<T,R> &a, Grid2<T,R> &b);	
	template <typename T, typename R> friend const T dot(const Grid2<T,R> &a, const Grid2<T,R> &b);	

public : 	
	void get_idx_normalized_coord(index_type &i, index_type &j, R& fx, R& fy, const Vector2<R>& _pos) const;

	void getIdx(index_type &i, index_type &j, const Vector2<R>& _pos) const;
	void getClosestIdx(index_type &i, index_type &j, const Vector2<R>& _pos) const;
	const T bilerp(index_type i, index_type j, R fx, R fy) const;

	inline void getIdx(Vector2<index_type>& idx, const Vector2<R>& _pos) const			{ getIdx(idx[0], idx[1], _pos);					}
	inline void getClosestIdx(Vector2<index_type>& idx, const Vector2<R>& _pos) const	{ getClosestIdx(idx[0], idx[1], _pos);			}
	inline const T bilerp(const Vector2<index_type>& idx, R fx, R fy) const				{ return bilerp(idx[0], idx[1], fx, fy);		}

	const bool exportBinary(std::ofstream&) const;
	const bool importBinary(std::ifstream&);

public : 	
	index_type dim[2];
	index_type nX, nY;
	index_type _size;

	T *data;

	R dx, dy;
	R dxInv, dyInv;
	R dxInvHalf, dyInvHalf;

	Vector2<R> offset, origin;
};

template<typename T, typename R> 
const bool Grid2<T,R>::exportBinary(std::ofstream& ofs) const {
	ofs.write((char*)(dim), sizeof(index_type) * 2);
	ofs.write((char*)(&nX), sizeof(index_type));
	ofs.write((char*)(&nY), sizeof(index_type));
	ofs.write((char*)(&_size), sizeof(index_type));

	ofs.write((char*)(data), sizeof(T) * _size);

	ofs.write((char*)(&dx), sizeof(R));		ofs.write((char*)(&dy), sizeof(R));
	ofs.write((char*)(&dxInv), sizeof(R));	ofs.write((char*)(&dyInv), sizeof(R));
	ofs.write((char*)(&dxInvHalf), sizeof(R));	ofs.write((char*)(&dyInvHalf), sizeof(R));

	ofs.write((char*)(&offset), sizeof(R)*2);
	ofs.write((char*)(&origin), sizeof(R)*2);

	return true;
}

template<typename T, typename R> 
const bool Grid2<T,R>::importBinary(std::ifstream& ifs) {
	ifs.read((char*)(dim), sizeof(index_type) * 2);
	ifs.read((char*)(&nX), sizeof(index_type));
	ifs.read((char*)(&nY), sizeof(index_type));
	ifs.read((char*)(&_size), sizeof(index_type));

	if(data) delete data; data = new T[_size];
	ifs.read((char*)(data), sizeof(T) * _size);

	ifs.read((char*)(&dx), sizeof(R));			ifs.read((char*)(&dy), sizeof(R));
	ifs.read((char*)(&dxInv), sizeof(R));		ifs.read((char*)(&dyInv), sizeof(R));
	ifs.read((char*)(&dxInvHalf), sizeof(R));	ifs.read((char*)(&dyInvHalf), sizeof(R));

	ifs.read((char*)(&offset), sizeof(R)*2);
	ifs.read((char*)(&origin), sizeof(R)*2);

	return true;
}

//////////////////////////////////////////////////////////////////////////
// Constructor/Destructor
template<typename T, typename R> 
Grid2<T,R>::Grid2() 
: nX(0), nY(0), _size(0), data(NULL), dx(0.0), dy(0.0), dxInv(0.0), dyInv(0.0), dxInvHalf(0.0), dyInvHalf(0.0), offset(0,0), origin(0,0) { dim[0] = dim[1] = 0; }

template<typename T, typename R> 
Grid2<T,R>::Grid2(const Grid2<T,R> &grid2)
: nX(0), nY(0), _size(0), data(NULL), dx(0.0), dy(0.0), dxInv(0.0), dyInv(0.0), dxInvHalf(0.0), dyInvHalf(0.0), offset(0,0), origin(0,0) {
	dim[0] = dim[1] = 0;
	resize(grid2);
	set(grid2);
}

template<typename T, typename R> 
Grid2<T,R>::Grid2(index_type _nX, index_type _nY, R _dx, R _dy, Vector2<R> _offset = Vector2<R>(0.0,0.0)) 
: nX(0), nY(0), _size(0), data(NULL), dx(0.0), dy(0.0), dxInv(0.0), dxInvHalf(0.0), dyInvHalf(0.0), dyInv(0.0), offset(0,0), origin(0,0) {
	dim[0] = dim[1] = 0;
	init(_nX, _nY, _dx, _dy, _offset); 
}

template<typename T, typename R> 
Grid2<T,R>::~Grid2() { 
	release(); 
}

template<typename T, typename R> void Grid2<T,R>::init(index_type _nX, index_type _nY, R _dx, R _dy, const Vector2<R> &_offset)	{
	release();
	dim[0] = _nX; dim[1] = _nY;
	nX = _nX;
	nY = _nY;
	_size = nX * nY;
	data = new T[_size];
	dx = _dx;	dy = _dy;
	
	if(IsAlmostZero(dx) || IsAlmostZero(dy)) { dxInv = 0; dyInv = 0; }
	else { dxInv = R(1.0)/dx; dyInv = R(1.0)/dy; dxInvHalf = R(0.5) * dxInv; dyInvHalf = R(0.5) * dyInv; }

	offset = _offset;
	origin = Vector2<R>(offset[0] * dx, offset[1] * dy);
}

template<typename T, typename R> template<typename S, typename R2> void Grid2<T,R>::resize(const Grid2<S,R2>& grid2) {
	init(grid2.nX, grid2.nY, grid2.dx, grid2.dy, Vector2<R>(grid2.offset));
}

template<typename T, typename R> void Grid2<T,R>::release()	{
	if(data) delete[] data; data = NULL;
	nX = nY = _size = 0;
}

template<typename T, typename R> void Grid2<T,R>::set_offset(const Vector2<R>& _offset) {
	offset = _offset;
	origin = Vector2<R>(offset[0] * dx, offset[1] * dy);
}

template<typename T, typename R> void Grid2<T,R>::set_origin(const Vector2<R>& _origin) {
	if(IsAlmostZero(dx) || IsAlmostZero(dy))
		basis::throwError("Should Initialize Grid2 first\n");

	origin = _origin;
	offset = Vector2<R>(origin[0] / dx, origin[1] / dy);
}

template<typename T, typename R> template<int DIM> const bool Grid2<T,R>::is_valid_range(index_type i) const {
	if(0 <= i && i < dim[DIM])	return true;
	else						return false;
}

template<typename T, typename R> template<int DIM> const bool Grid2<T,R>::is_valid_range_for_interp(index_type i) const {
	if(0 <= i && i < dim[DIM]-1)	return true;
	else							return false;
}

template<typename T, typename R> const bool Grid2<T,R>::is_empty() const {
	if(nX == 0 || nY == 0)
		return true;
	return false;
}

//////////////////////////////////////////////////////////////////////////
// Index operators
template<typename T, typename R> inline const T Grid2<T,R>::operator[](index_type i) const				 { __ASSERT(i<_size); return data[i];			}
template<typename T, typename R> inline const T Grid2<T,R>::operator()(index_type i) const				 { __ASSERT(i<_size); return data[i];			}
template<typename T, typename R> inline const T Grid2<T,R>::operator()(index_type i, index_type j) const { __ASSERT(i<nX && j<nY); return data[i*nY+j];	}

template<typename T, typename R> inline T &Grid2<T,R>::operator[](index_type i)					{ __ASSERT(i<_size); return data[i];			}
template<typename T, typename R> inline T &Grid2<T,R>::operator()(index_type i)					{ __ASSERT(i<_size); return data[i];			}
template<typename T, typename R> inline T &Grid2<T,R>::operator()(index_type i, index_type j)	{ __ASSERT(i<nX && j<nY); return data[i*nY+j];	}

template<typename T, typename R> inline const T Grid2<T,R>::operator()(const Vector2<index_type>& i) const { __ASSERT(i[0]<nX && i[1]<nY); return data[i[0]*nY+i[1]];	}
template<typename T, typename R> inline T &Grid2<T,R>::operator()(const Vector2<index_type>& i)			   { __ASSERT(i[0]<nX && i[1]<nY); return data[i[0]*nY+i[1]];	}

template<typename T, typename R> inline const T Grid2<T,R>::get(index_type i, index_type j) const   { __ASSERT(i<nX && j<nY);		return data[i*nY+j];		}
template<typename T, typename R> inline const T Grid2<T,R>::get(const Vector2<index_type>& i) const { __ASSERT(i[0]<nX && i[1]<nY); return data[i[0]*nY+i[1]];	}

//////////////////////////////////////////////////////////////////////////
// Getter
template<typename T, typename R> inline const T Grid2<T,R>::get(const Vector2<R> &_pos) const {
	return (*this)(_pos);
}

template<typename T, typename R> inline const T Grid2<T,R>::operator()(const Vector2<R> &_pos) const { 
	index_type i,j;
	getIdx(i,j,_pos);
	const Vector2<R> _corner = pos(i,j);
	const R fx = (_pos[0] - _corner[0]) * dxInv;
	const R fy = (_pos[1] - _corner[1]) * dyInv;
	return bilerp(i,j,fx,fy);	   
}

template<typename T, typename R> inline const T Grid2<T,R>::get_nearest(const Vector2<R>& _pos) const {
	index_type i,j;
	getIdx(i,j,_pos);
	Vector2<R> _corner = pos(i,j);
	const R fx = (_pos[0] - _corner[0]) * dxInv;
	const R fy = (_pos[1] - _corner[1]) * dyInv;

	index_type i_icr = (fx < 0.5) ? 0 : 1;
	index_type j_icr = (fy < 0.5) ? 0 : 1;

//	i += i_icr; j += j_icr;

//	if(i < 0) i = 0; else if(i > nX-1) i = nX-1;
//	if(j < 0) j = 0; else if(j > nY-1) j = nY-1;

	return (*this)(i+i_icr,j+j_icr);
}

template<typename T, typename R> inline const Vector2<R> Grid2<T,R>::pos(index_type i, index_type j) const	 { return Vector2<R>(i*dx+origin[0], j*dy+origin[1]);		}
template<typename T, typename R> inline const Vector2<R> Grid2<T,R>::pos(const index_type i[2]) const { return Vector2<R>(i[0]*dx+origin[0], i[1]*dy+origin[1]);	}
template<typename T, typename R> inline const Vector2<R> Grid2<T,R>::pos(const Vector2<index_type>& i) const { return Vector2<R>(i[0]*dx+origin[0], i[1]*dy+origin[1]);	}

template<typename T, typename R> inline void Grid2<T,R>::getMaxAndMin(T &_max, T &_min, const Vector2<R> &_pos) const {
	index_type i,j;
	getIdx(i,j,_pos);
	_max = MAX( MAX((*this)(i,j),(*this)(i+1,j)), MAX((*this)(i,j+1),(*this)(i+1,j+1)) );
	_min = MIN( MIN((*this)(i,j),(*this)(i+1,j)), MIN((*this)(i,j+1),(*this)(i+1,j+1)) );
}

//////////////////////////////////////////////////////////////////////////
// Assignment and set
template<typename T, typename R> void Grid2<T,R>::zero()					{ VectorOperation<Grid2<T,R>>().zero(*this); }

template<typename T, typename R> template<typename T1, typename R1> void Grid2<T,R>::set(const Grid2<T1,R1> &grid) {
	if(nX != grid.nX || nY != grid.nY)
		throwError("Error in Grid2<T,R>::void set(const Grid2<T1,R1> &grid)");

	for(index_type i=0;i<nX;++i)
		for(index_type j=0;j<nY;++j)
			(*this)(i,j) = grid(i,j);
}
template<typename T, typename R> void Grid2<T,R>::set(const Grid2& grid)	{ VectorOperation<Grid2<T,R>>().set(*this, grid);	}
template<typename T, typename R> void Grid2<T,R>::set(const T val)			{ VectorOperation<Grid2<T,R>>().set(*this, val);	}

template<typename T, typename R> void Grid2<T,R>::set_boundary(const T val) {
	for(index_type i=0;i<nX;i++) (*this)(i,0) = (*this)(i,nY-1) = val;
	for(index_type j=0;j<nY;j++) (*this)(0,j) = (*this)(nX-1,j) = val;
}

template<typename T, typename R> const Grid2<T,R> &Grid2<T,R>::operator=(const Grid2<T,R> &grid) { set(grid); return (*this); }

//////////////////////////////////////////////////////////////////////////
// Manipulation	
template <typename T, typename R> void Grid2<T,R>::setNeumannBoundary() {
	GridBoundaryCondition_Neumann<T,R> neumann; neumann.enforce(*this);	
}

template <typename T, typename R> const T Grid2<T,R>::infNorm() const {
	VectorFunctor_InfNorm<Grid2> func(*this);
	func(blocked_range<index_type>(0,_size));
	return func.getValue();	
}

template <typename T, typename R> const T Grid2<T,R>::l2Norm() const {
	VectorFunctor_L2Norm<Grid2> func(*this);
	func(blocked_range<index_type>(0,_size));
	return func.getValue();		
}

template <typename T, typename R> const T Grid2<T,R>::sum() const {
	VectorFunctor_Sum<Grid2> func(*this);
	func(blocked_range<index_type>(0,_size));
	return func.getValue();		
}

//////////////////////////////////////////////////////////////////////////
// In place arithmetic
template<typename T, typename R> const Grid2<T,R> &Grid2<T,R>::operator+=(const Grid2 &g) { VectorFunctor_Icr<Grid2>(*this, g)(blocked_range<index_type>(0,_size)); return (*this);	}
template<typename T, typename R> const Grid2<T,R> &Grid2<T,R>::operator-=(const Grid2 &g) { VectorFunctor_Dcr<Grid2>(*this, g)(blocked_range<index_type>(0,_size));	return (*this);	}
template<typename T, typename R> const Grid2<T,R> &Grid2<T,R>::operator*=(const R c)	  { VectorFunctor_Scale<Grid2,R>(*this, c)(blocked_range<index_type>(0,_size));return (*this); }	

template<typename T, typename R> inline const geometryspace::AABB2D<R> Grid2<T,R>::getAABB() const {
	return geometryspace::AABB2D<R>(origin, origin + Vector2<R>(dx * (nX-1), dy * (nY-1)));
}

//////////////////////////////////////////////////////////////////////////
// glDraw
template<typename T, typename R> void Grid2<T,R>::glDrawGrid(const RGBA &color, GLfloat width) const {
	glLineWidth(width);
	glColor4f( color[0], color[1], color[2], color[3] );

	R Lx = dx * (nX-1), Ly = dy * (nY-1);

	for(index_type i=0;i<nX;++i)	{
		glBegin( GL_LINES );
		glVertex3f( (GLfloat)(i*dx + origin[0]), (GLfloat)origin[1]			, 0.0f);
		glVertex3f( (GLfloat)(i*dx + origin[0]), (GLfloat)(origin[1] + Ly)	, 0.0f);
		glEnd();		
	}
	for(index_type i=0;i<nY;++i)	{
		glBegin( GL_LINES );
		glVertex3f( (GLfloat)origin[0]		 , (GLfloat)(i*dy + origin[1]) , 0.0f);
		glVertex3f( (GLfloat)(origin[0] + Lx), (GLfloat)(i*dy + origin[1]) , 0.0f);
		glEnd();		
	}
}

template<typename T, typename R> void Grid2<T,R>::save_as_matlab(const char* const filename, const char* const grid_name) const {
	std::ofstream ofs(filename);
	if(!ofs) { 
		throwError("Error in void Grid2<T,R>::save_as_matlab(const char* const filename, const char* const grid_name)");
		return;
	}
	ofs << grid_name << " = ones(" << nX-2 << "," << nY-2 << ");" << std::endl;
	for(index_type i=1;i<nX-1;++i) {
		for(index_type j=1;j<nY-1;++j) {
			ofs << grid_name << "(" << i << "," << j << ") = " << (*this)(i,j) << ";" << std::endl;
		}
	}
}

//////////////////////////////////////////////////////////////////////////
// Private Functions
template<typename T, typename R> void Grid2<T,R>::getIdx(index_type &i, index_type &j, const Vector2<R>& _pos) const {
	i = (index_type)((_pos[0] - origin[0]) * dxInv);
	j = (index_type)((_pos[1] - origin[1]) * dyInv);

	if(i == nX-1) --i;	if(j == nY-1) --j;
}

template<typename T, typename R> void Grid2<T,R>::get_idx_normalized_coord(index_type &i, index_type &j, R& fx, R& fy, const Vector2<R>& _pos) const {
	getIdx(i,j,_pos);
	Vector2<R> _corner = pos(i,j);
	fx = (_pos[0] - _corner[0]) * dxInv;
	fy = (_pos[1] - _corner[1]) * dyInv;
}

template<typename T, typename R> void Grid2<T,R>::getClosestIdx(index_type &i, index_type &j, const Vector2<R>& _pos) const {
	R ir = ((_pos[0] - origin[0]) * dxInv);
	R jr = ((_pos[1] - origin[1]) * dyInv);

	i = (index_type)ir; j = (index_type)jr;

	if(ir - i > 0.5f) ++i;
	if(jr - j > 0.5f) ++j;

	if(i == nX-1) --i;	else if(i == 0) ++i;
	if(j == nY-1) --j;	else if(j == 0)	++j;
}

template<typename T, typename R> const T Grid2<T,R>::bilerp(index_type i, index_type j, R fx, R fy) const {
	return (1-fx)*((1-fy)*(*this)(i,j)+fy*(*this)(i,j+1))+fx*((1-fy)*(*this)(i+1,j)+fy*(*this)(i+1,j+1)); 
}

template <typename T, typename R> inline void swap(Grid2<T,R> &a, Grid2<T,R> &b)	{
	T *tmp = a.data;	a.data = b.data;	b.data = tmp;	
}

template <typename T, typename R> const T dot(const Grid2<T,R> &a, const Grid2<T,R> &b)	{
	ArrayDotFunctor<T> dotfunc(a.data, b.data);
	dotfunc(blocked_range<typename Grid2<T,R>::index_type>(0, a._size));
	return dotfunc.getValue();
}

#undef __ASSERT

};	// namespace arithmetic

#endif