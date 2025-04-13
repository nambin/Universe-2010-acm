#ifndef _ARITHMETIC_MATRIXN_H_
#define _ARITHMETIC_MATRIXN_H_

#include "../../../Basis/Basis/BasicLibrary.h"
#include "../Vector/VectorN.h"
#include "../MathLibrary/MatrixVectorOperation.h"
#include "MatrixOperation.h"

namespace arithmetic {

#ifdef _DEBUG
#define __ASSERT(c) assert(c)
#else
#define __ASSERT(c) {}
#endif

//////////////////////////////////////////////////////////////////////////
//	MatrixN Class Definition.
//	NxM Matrix Class.
//////////////////////////////////////////////////////////////////////////
template <typename T> class MatrixN {
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
	// Constructor/Destructor		
	MatrixN(index_type n = 0, index_type m = 0, const T val = T(0) );

	MatrixN(const MatrixN<T> &src );
	MatrixN(index_type n, index_type m, const T **a);
	MatrixN(index_type n, index_type m, const VectorN<T>* v );		

	~MatrixN();

	void init(index_type n = 0, index_type m = 0);

	//////////////////////////////////////////////////////////////////////////
	// Index operators
	inline const T  get(index_type i, index_type j) const;		

	inline T&		operator()(index_type i, index_type j);		
	inline const T  operator()(index_type i, index_type j) const;

	//////////////////////////////////////////////////////////////////////////
	// Getter
	void getRow(index_type row, VectorN<T> &v) const;
	void getCol(index_type col, VectorN<T> &v) const;

	//////////////////////////////////////////////////////////////////////////
	// Size
	inline const index_type dimRow() const	{ return _dimRow; }
	inline const index_type dimCol() const	{ return _dimCol; } 
	inline const index_type nRows() const	{ return _dimRow; }
	inline const index_type nCols() const	{ return _dimCol; } 

	//////////////////////////////////////////////////////////////////////////
	// Resize Matrix
	void resize( index_type n, index_type m );

	//////////////////////////////////////////////////////////////////////////
	// Assignment and set
	void zero();
	void identity();

	inline void set( index_type i, index_type j, const T d) ;
	void setRow( index_type row, const VectorN<T> &v);
	void setCol( index_type col, const VectorN<T> &v);

	void set( const MatrixN<T> &src );
	void set( const T **src ); 
	void set( const VectorN<T>* v);
	void set( const T val);

	MatrixN &operator=( const MatrixN &src );
	MatrixN &operator=( const T **src );
	MatrixN &operator=( const VectorN<T>* v);
	MatrixN &operator=( const Matrix_neg_M<MatrixN<T>> & op);
	MatrixN &operator=( const Matrix_trans_M<MatrixN<T>> & op);
	MatrixN &operator=( const Matrix_M_add_M<MatrixN<T>> & op);
	MatrixN &operator=( const Matrix_M_dif_M<MatrixN<T>> & op);
	MatrixN &operator=( const Matrix_M_mul_S<MatrixN<T>,T> & op);
	MatrixN &operator=( const Matrix_M_mul_M<MatrixN<T>> & op);

	//////////////////////////////////////////////////////////////////////////
	// Manipulation
	const T frobeniusNorm() const;
	MatrixN &transpose();

	//////////////////////////////////////////////////////////////////////////
	// Get Pointer
	inline const T* const getptr() const { return data; }

	//////////////////////////////////////////////////////////////////////////
	// Query
	const bool isZero() const ;
	const bool isIdentity() const ;
	const bool isSquareMatrix() const	{ return (_dimRow == _dimCol);	}	
	const bool isSymmetricMatrix() const;

	//////////////////////////////////////////////////////////////////////////
	// In place arithmetic
	MatrixN &operator+=(T d);	
	MatrixN &operator-=(T d);	
	MatrixN &operator*=(T d);	
	MatrixN &operator/=(T d);	

	//////////////////////////////////////////////////////////////////////////
	// In place Arithmetic Operator
	MatrixN &operator+=(const MatrixN& a);
	MatrixN &operator-=(const MatrixN& a);

	//////////////////////////////////////////////////////////////////////////
	// Row operation
	MatrixN &scaleRow( index_type row, const T c );
	MatrixN &addRow( index_type row0, index_type row1, const T c );
	MatrixN &interchangeRow( index_type row0, index_type row1 );

	//////////////////////////////////////////////////////////////////////////
	// IO
	template <typename T> friend std::ostream& operator<<(std::ostream &os, const MatrixN<T> &m);	

	//////////////////////////////////////////////////////////////////////////
	// Friends Operators
	template <typename T> friend void neg( MatrixN<T> &ret, const MatrixN<T> &a);
	template <typename T> friend void add( MatrixN<T> &ret, const MatrixN<T> &a, const MatrixN<T> &b);
	template <typename T> friend void dif( MatrixN<T> &ret, const MatrixN<T> &a, const MatrixN<T> &b);

	template <typename T> friend void mul( MatrixN<T> &ret, const MatrixN<T> &m, const T &c);
	template <typename T> friend void mul( MatrixN<T> &ret, const T &c, const MatrixN<T> &m );

	template <typename T> friend void div( MatrixN<T> &ret, const MatrixN<T> &m, const T &c);

	template <typename T> friend const T frobeniusNorm( const MatrixN<T> &a );

	template <typename T> friend void mul( VectorN<T> &ret, const MatrixN<T> &m, const VectorN<T> &v);
	template <typename T> friend void mul( MatrixN<T> &ret, const MatrixN<T> &a, const MatrixN<T> &b);
	template <typename T> friend void transpose( MatrixN<T> &trans, const MatrixN<T> &a);

	template <typename T> friend const Matrix_neg_M<MatrixN<T>> operator-(const MatrixN<T> &a);
	template <typename T> friend const Matrix_trans_M<MatrixN<T>> transpose(const MatrixN<T> &a);
	template <typename T> friend const Matrix_M_add_M<MatrixN<T>> operator+(const MatrixN<T> &a, const MatrixN<T> &b);
	template <typename T> friend const Matrix_M_dif_M<MatrixN<T>> operator-(const MatrixN<T> &a, const MatrixN<T> &b);
	template <typename T> friend const Matrix_M_mul_S<MatrixN<T>,T> operator*(const MatrixN<T> &a, const T &c);
	template <typename T> friend const Matrix_M_mul_S<MatrixN<T>,T> operator*(const T &c, const MatrixN<T> &a);
	template <typename T> friend const Matrix_M_mul_S<MatrixN<T>,T> operator/(const MatrixN<T> &a, const T &c);
	template <typename T> friend const Matrix_M_mul_M<MatrixN<T>> operator*(const MatrixN<T> &a, const MatrixN<T> &b);

	template <typename T> friend const Vector_M_mul_V<MatrixN<T>, VectorN<T>> operator*(const MatrixN<T> &m, const VectorN<T> &v);
	template <typename T> friend const Vector_S_mul_M_mul_V<T, MatrixN<T>, VectorN<T>> operator*(const T s, const Vector_M_mul_V<MatrixN<T>, VectorN<T>> &Mv);
	template <typename T> friend const Vector_V_sum_M_mul_V<MatrixN<T>, VectorN<T>> operator+(const VectorN<T> &v, const Vector_M_mul_V<MatrixN<T>, VectorN<T>> &Mv);
	template <typename T> friend const Vector_V_sub_M_mul_V<MatrixN<T>, VectorN<T>> operator-(const VectorN<T> &v, const Vector_M_mul_V<MatrixN<T>, VectorN<T>> &Mv);
	template <typename T> friend const Vector_V_sum_S_mul_M_mul_V<T, MatrixN<T>, VectorN<T>> operator+(const VectorN<T> &v, const Vector_S_mul_M_mul_V<T, MatrixN<T>, VectorN<T>> &sMv);
	template <typename T> friend const Vector_V_sub_S_mul_M_mul_V<T, MatrixN<T>, VectorN<T>> operator-(const VectorN<T> &v, const Vector_S_mul_M_mul_V<T, MatrixN<T>, VectorN<T>> &sMv);

private : 	
	inline const index_type idx(index_type row, index_type col) const { return (row * _dimCol + col);	}

private :
	index_type	_dimRow, _dimCol;
	T *			data;	
};

//////////////////////////////////////////////////////////////////////////
// Constructor/Destructor
template <typename T> MatrixN<T>::MatrixN(index_type n, index_type m, const T val ) : data(NULL) {
	__ASSERT(n>=0 && m>=0);
	if(n == 0 || m == 0) {
		_dimRow = _dimCol = 0;
		data = NULL;
		return;
	}

	init(n, m);
	set(val);
}

template <typename T> MatrixN<T>::MatrixN( const MatrixN<T> &src ): data(NULL) {
	init(src.dimRow(), src.dimCol());
	set(src);
}

template <typename T> MatrixN<T>::MatrixN( index_type n, index_type m, const T **a ): data(NULL) {
	init(n, m);
	set(a);
}

template <typename T> MatrixN<T>::MatrixN( index_type n, index_type m, const VectorN<T>* v ): data(NULL) {
	__ASSERT(m == v[0].dim());
	init(n, m);
	set(v);
}

template <typename T> MatrixN<T>::~MatrixN() {
	if(data) delete[] data;
}

template <typename T> void MatrixN<T>::init(index_type n, index_type m) {
	__ASSERT(n>=0 && m>=0);

	if(data) delete[] data;

	if(n == 0 || m == 0) {
		_dimRow = _dimCol = 0;
		data = NULL;
	}
	else {
		_dimRow = n;	
		_dimCol = m;
		data = new T[_dimRow * _dimCol];
		zero();
	}
}

//////////////////////////////////////////////////////////////////////////
// Index operators
template <typename T> inline const T MatrixN<T>::get(index_type i, index_type j) const {
	assert( 0 <= i && 0 <= j && i < _dimRow && j < _dimCol);
	return data[idx(i,j)];
}

template <typename T> inline T&	MatrixN<T>::operator()(index_type i, index_type j) {
	__ASSERT( 0 <= i && 0 <= j && i < _dimRow && j < _dimCol);
	return data[idx(i,j)];
}

template <typename T> inline const T MatrixN<T>::operator()(index_type i, index_type j) const {
	__ASSERT( 0 <= i && 0 <= j && i < _dimRow && j < _dimCol);
	return data[idx(i,j)];
}

//////////////////////////////////////////////////////////////////////////
// Getter
template <typename T> void MatrixN<T>::getRow(index_type row, VectorN<T> &v) const {
	__ASSERT(v.dim() == _dimCol);

	for(index_type i=0;i<_dimCol;++i)
		v[i] = (*this)(row, i);
}

template <typename T> void MatrixN<T>::getCol(index_type col, VectorN<T> &v) const {
	__ASSERT(v.dim() == _dimRow);

	for(index_type i=0;i<_dimRow;++i)
		v[i] = (*this)(i, col);
}

//////////////////////////////////////////////////////////////////////////
// Resize Matrix
template <typename T> void MatrixN<T>::resize(index_type n, index_type m)	{	
	init(n,m);	
}

//////////////////////////////////////////////////////////////////////////
// Assignment and set// Simple Setter
template <typename T> void MatrixN<T>::zero()	{
	for(index_type i=0;i<dimRow();++i) 
		for(index_type j=0;j<dimCol();++j)
			(*this)(i,j) = 0;
}

template <typename T> void MatrixN<T>::identity()	{
	__ASSERT(_dimRow == _dimCol);

	zero();
	for(index_type i=0;i<dimRow();++i)
		(*this)(i,i) = T(1);
}


template <typename T> inline void MatrixN<T>::set( index_type i, index_type j, T d)	{
	__ASSERT(  0 <= i && 0 <= j && i < _dimRow && j < _dimCol );
	(*this)(i,j) = d;
}

template <typename T> void MatrixN<T>::setRow(index_type row, const VectorN<T> &v) {
	__ASSERT( 0 <= row && row < _dimRow && _dimCol == v.size());

	for(index_type i=0;i<_dimCol;++i)
		(*this)(row, i) = v(i);
}

template <typename T> void MatrixN<T>::setCol(index_type col, const VectorN<T> &v) {
	__ASSERT( 0 <= col && col < _dimCol && _dimRow == v.size());

	for(index_type i=0;i<_dimRow;++i)
		(*this)(i, col) = v(i);
}

template <typename T> void MatrixN<T>::set( const MatrixN<T> &src )	{
	__ASSERT( _dimRow == src._dimRow && _dimCol == src._dimCol );

	for(index_type i=0;i<dimRow();++i) 
		for(index_type j=0;j<dimCol();++j)
			(*this)(i,j) = src(i,j);
}

template <typename T> void MatrixN<T>::set( const T **a )	{
	for(index_type i=0;i<dimRow();++i) 
		for(index_type j=0;j<dimCol();++j)
			(*this)(i,j) = a[i][j];
}

template <typename T> void MatrixN<T>::set( const VectorN<T>* v )	{
	for(index_type i=0;i<dimRow();++i) 
		for(index_type j=0;j<dimCol();++j)
			(*this)(i,j) = v[i][j];
}

template <typename T> void MatrixN<T>::set( const T val )	{
	for(index_type i=0;i<dimRow();++i) 
		for(index_type j=0;j<dimCol();++j)
			(*this)(i,j) = val;
}

template <typename T> MatrixN<T> &MatrixN<T>::operator=( const MatrixN &src )	{	
	set(src);
	return (*this);
}

template <typename T> MatrixN<T> &MatrixN<T>::operator=( const T **a )	{	
	set(a);
	return (*this);
}

template <typename T> MatrixN<T> &MatrixN<T>::operator=( const VectorN<T>* v )	{	
	set(v);
	return (*this);
}

template <typename T> MatrixN<T> &MatrixN<T>::operator=( const Matrix_neg_M<MatrixN<T>> & op)	{
	arithmetic::neg(*this, op.M);
	return (*this);
}

template <typename T> MatrixN<T> &MatrixN<T>::operator=( const Matrix_trans_M<MatrixN<T>> & op)	{
	arithmetic::transpose(*this, op.M);
	return (*this);
}

template <typename T> MatrixN<T> &MatrixN<T>::operator=( const Matrix_M_add_M<MatrixN<T>> & op)	{
	arithmetic::add(*this, op.M1, op.M2);
	return (*this);
}

template <typename T> MatrixN<T> &MatrixN<T>::operator=( const Matrix_M_dif_M<MatrixN<T>> & op)	{
	arithmetic::dif(*this, op.M1, op.M2);
	return (*this);
}

template <typename T> MatrixN<T> &MatrixN<T>::operator=( const Matrix_M_mul_S<MatrixN<T>,T> & op)	{
	arithmetic::mul(*this, op.M, op.s);
	return (*this);
}

template <typename T> MatrixN<T> &MatrixN<T>::operator=( const Matrix_M_mul_M<MatrixN<T>> & op)	{
	arithmetic::mul(*this, op.M1, op.M2);
	return (*this);
}

//////////////////////////////////////////////////////////////////////////
// Manipulation
template <typename T> const T MatrixN<T>::frobeniusNorm() const {
	return arithmetic::frobeniusNorm(*this);
}

template <typename T> MatrixN<T>& MatrixN<T>::transpose() {
	__ASSERT(_dimRow == _dimCol);

	T tmp;
	for(index_type i=0; i<_dimRow ; ++i) {
		for(index_type j=0; j<i ; ++j) {
			tmp = (*this)(i,j);
			(*this)(i,j) = (*this)(j,i);
			(*this)(j,i) = tmp;
		}
	}

	return (*this);
}

//////////////////////////////////////////////////////////////////////////
// Query
template <typename T> const bool MatrixN<T>::isZero() const	{
	for(index_type i=0;i<dimRow();++i) 
		for(index_type j=0;j<dimCol();++j)
			if(!IsAlmostZero((*this)(i,j)))
				return false;
	return true;
}

template <typename T> const bool MatrixN<T>::isIdentity() const	{
	__ASSERT(_dimRow == _dimCol);

	for(index_type i=0;i<dimRow();++i) {
		for(index_type j=0;j<dimCol();++j)	{
			if(i==j && !IsAlmostZero((*this)(i,j)-1))
				return false;
			if(i!=j && !IsAlmostZero((*this)(i,j)))
				return false;
		}
	}
	return true;	
}

template <typename T> const bool MatrixN<T>::isSymmetricMatrix() const {
	if(!isSquareMatrix())
		return false;

	for(index_type i=0;i<dimRow();++i) {
		for(index_type j=i+1;j<dimCol();++j) {
			if(!IsAlmostZero((*this)(i,j) - (*this)(j,i)))
				return false;
		}
	}
	return true;
}

//////////////////////////////////////////////////////////////////////////
// In place arithmetic
template <typename T> MatrixN<T> &MatrixN<T>::operator+=(T d) {
	for(index_type i=0;i<dimRow();++i) 
		for(index_type j=0;j<dimCol();++j)	
			(*this)(i,j) += d;
	return (*this);
}

template <typename T> MatrixN<T> &MatrixN<T>::operator-=(T d) {
	return (*this) += (-d);
}

template <typename T> MatrixN<T> &MatrixN<T>::operator*=(T d) {
	for(index_type i=0;i<dimRow();++i) 
		for(index_type j=0;j<dimCol();++j)	
			(*this)(i,j) *= d;
	return (*this);
}

template <typename T> MatrixN<T> &MatrixN<T>::operator/=(T d) {
	return (*this) *= (1.0f / d);
}

//////////////////////////////////////////////////////////////////////////
// In place Arithmetic Operator
template <typename T> MatrixN<T> &MatrixN<T>::operator+=( const MatrixN& a)	{
	__ASSERT( _dimRow == a._dimRow && _dimCol == a._dimCol );

	for(index_type i=0;i<dimRow();++i) 
		for(index_type j=0;j<dimCol();++j)
			(*this)(i,j) += a(i,j);	

	return (*this);
}

template <typename T> MatrixN<T> &MatrixN<T>::operator-=( const MatrixN& a)	{
	__ASSERT( ( _dimRow == a._dimRow && _dimCol == a._dimCol ) );

	for(index_type i=0;i<dimRow();++i) 
		for(index_type j=0;j<dimCol();++j)
			(*this)(i,j) -= a(i,j);	

	return (*this);
}

//////////////////////////////////////////////////////////////////////////
// Row operation
template <typename T> MatrixN<T> &MatrixN<T>::scaleRow( index_type row, const T c )	{
	for(index_type i=0;i<dimCol();++i)
		(*this)(row,i) *= c;
	return (*this);
}

template <typename T> MatrixN<T> &MatrixN<T>::addRow( index_type row0, index_type row1, const T c )	{
	for(index_type i=0;i<dimCol();++i)
		(*this)(row0,i) += (*this)(row1,i) * c;
	return (*this);
}

template <typename T> MatrixN<T> &MatrixN<T>::interchangeRow( index_type row0, index_type row1 )	{
	T tmp;
	for(index_type i=0;i<dimCol();++i) {
		tmp = (*this)(row0, i);
		(*this)(row0, i) = (*this)(row1, i);
		(*this)(row1, i) = tmp;
	}
	return (*this);
}

//////////////////////////////////////////////////////////////////////////
// IO
template <typename T> std::ostream& operator<<(std::ostream &os, const MatrixN<T> &m) {
	for(MatrixN<T>::index_type i=0;i<m.dimRow();++i) {
		for(MatrixN<T>::index_type j=0;j<m.dimCol();++j) {
			os << std::setw(8) << m(i,j) << " ";			
		}
		os << std::endl;		
	}
	os << std::endl;
	return os;
}


//////////////////////////////////////////////////////////////////////////
// Friends Functions
template <typename T> void neg( MatrixN<T> &ret, const MatrixN<T> &a)	{
	__ASSERT( ( ret.dimRow() == a.dimRow() && ret.dimCol() == a.dimCol() ) );

	for(index_type i=0; i<ret.dimRow(); ++i)	
		for(index_type j=0; j<ret.dimCol(); ++j)	
			ret(i,j) = -a(i,j);
}

template <typename T> void add( MatrixN<T> &ret, const MatrixN<T> &a, const MatrixN<T> &b)	{
	__ASSERT( ( ret.dimRow() == a.dimRow() && ret.dimCol() == a.dimCol() )
		   && ( ret.dimRow() == b.dimRow() && ret.dimCol() == b.dimCol() ));

	for(index_type i=0; i<ret.dimRow(); ++i)	
		for(index_type j=0; j<ret.dimCol(); ++j)	
			ret(i,j) = a(i,j) + b(i,j);			
}

template <typename T> void dif( MatrixN<T> &ret, const MatrixN<T> &a, const MatrixN<T> &b)	{
	__ASSERT( ( ret.dimRow() == a.dimRow() && ret.dimCol() == a.dimCol() )
		   && ( ret.dimRow() == b.dimRow() && ret.dimCol() == b.dimCol() ));

	for(index_type i=0; i<ret.dimRow(); ++i)	
		for(index_type j=0; j<ret.dimCol(); ++j)	
			ret(i,j) = a(i,j) - b(i,j);			
}

template <typename T> void mul( MatrixN<T> &ret, const MatrixN<T> &m, const T &c)	{
	__ASSERT( ( ret.dimRow() == m.dimRow() && ret.dimCol() == m.dimCol() ) );

	for(index_type i=0; i<ret.dimRow(); ++i)	
		for(index_type j=0; j<ret.dimCol(); ++j)	
			ret(i,j) = m(i,j) * c;	
}

template <typename T> void mul( MatrixN<T> &ret, const T &c, const MatrixN<T> &m )	{
	mul(ret,m,c);
}

template <typename T> void div( MatrixN<T> &ret, const MatrixN<T> &m, const T &c)	{
	mul(ret, m, T(1.0f) / c);	
}

template <typename T> const T frobeniusNorm( const MatrixN<T> &a ) {
	T sqrmag = (T)0;
	for(MatrixN<T>::index_type i=0; i<a.dimRow() ; ++i)			
		for(MatrixN<T>::index_type j=0; j<a.dimCol() ; ++j) 
			sqrmag += a(i,j) * a(i,j);
	return (T)sqrt(sqrmag);
}

template <typename T> void mul( VectorN<T> &ret, const MatrixN<T> &m, const VectorN<T> &v)	{
	__ASSERT( ret.dim() == m.dimRow() && v.dim() == m.dimCol() );

	if(&ret == &v) {
		VectorN<T> _ret(ret.dim());

		mul(_ret, m, v);
		ret = _ret;

		return;
	}

	ret.zero();
	for(MatrixN<T>::index_type i=0; i<m.dimRow() ; ++i)			
		for(MatrixN<T>::index_type j=0; j<m.dimCol() ; ++j) 
			ret(i) += m(i,j) * v(j);
}

template <typename T> void mul( MatrixN<T> &ret, const MatrixN<T> &a, const MatrixN<T> &b)	{
	__ASSERT( ret.dimRow() == a.dimRow() && ret.dimCol() == b.dimCol() && a.dimCol() == b.dimRow() );

	if(&ret == &a || &ret == &b) {
		MatrixN<T> _ret(ret.dimRow(), ret.dimCol());
		
		mul(_ret, a, b);
		ret = _ret;

		return;
	}

	ret.zero();
	for(MatrixN<T>::index_type i=0; i<ret.dimRow() ; ++i)			
		for(MatrixN<T>::index_type j=0; j<ret.dimCol() ; ++j) 
			for(MatrixN<T>::index_type k=0; k<a.dimCol(); ++k)	
				ret(i,j) += (a(i,k) * b(k,j));			
}

template <typename T> void transpose( MatrixN<T> &trans, const MatrixN<T> &a)	{
	__ASSERT( (trans.dimRow() == a.dimCol() && trans.dimCol() == a.dimRow()) );

	if(&trans == &a) {		
		trans.transpose();
		return;
	}

	for(MatrixN<T>::index_type i=0; i<trans.dimRow() ; ++i)			
		for(MatrixN<T>::index_type j=0; j<trans.dimCol() ; ++j) 
			trans(i,j) = a(j,i);
}

//////////////////////////////////////////////////////////////////////////
// Matrix Operations
template <typename T> const Matrix_neg_M<MatrixN<T>> operator-(const MatrixN<T> &a) {
	return Matrix_neg_M<MatrixN<T>>(a);
}

template <typename T> const Matrix_trans_M<MatrixN<T>> transpose(const MatrixN<T> &a) {
	return Matrix_trans_M<MatrixN<T>>(a);
}

template <typename T> const Matrix_M_add_M<MatrixN<T>> operator+(const MatrixN<T> &a, const MatrixN<T> &b) {
	return Matrix_M_add_M<MatrixN<T>>(a,b);
}

template <typename T> const Matrix_M_dif_M<MatrixN<T>> operator-(const MatrixN<T> &a, const MatrixN<T> &b) {
	return Matrix_M_dif_M<MatrixN<T>>(a,b);
}

template <typename T> const Matrix_M_mul_S<MatrixN<T>, T> operator*(const MatrixN<T> &a, const T &c) {
	return Matrix_M_mul_S<MatrixN<T>, T>(a,c);
}

template <typename T> const Matrix_M_mul_S<MatrixN<T>, T> operator*(const T &c, const MatrixN<T> &a) {
	return Matrix_M_mul_S<MatrixN<T>, T>(a,c);
}

template <typename T> const Matrix_M_mul_S<MatrixN<T>, T> operator/(const MatrixN<T> &a, const T &c) {
	return Matrix_M_mul_S<MatrixN<T>, T>(a,T(1.0f)/c);
}

template <typename T> const Matrix_M_mul_M<MatrixN<T>> operator*(const MatrixN<T> &a, const MatrixN<T> &b) {
	return Matrix_M_mul_M<MatrixN<T>>(a,b);
}

//////////////////////////////////////////////////////////////////////////
// Matrix/Vector Operations
template <typename T> const Vector_M_mul_V<MatrixN<T>, VectorN<T>> operator*(const MatrixN<T> &m, const VectorN<T> &v) {
	return Vector_M_mul_V<MatrixN<T>, VectorN<T>>(m,v);
}

template <typename T> const Vector_S_mul_M_mul_V<T, MatrixN<T>, VectorN<T>> operator*(const T s, const Vector_M_mul_V<MatrixN<T>, VectorN<T>> &Mv) {
	return Vector_S_mul_M_mul_V<T, MatrixN<T>, VectorN<T>>(s, Mv.M, Mv.a);
}

template <typename T> const Vector_V_sum_M_mul_V<MatrixN<T>, VectorN<T>> operator+(const VectorN<T> &v, const Vector_M_mul_V<MatrixN<T>, VectorN<T>> &Mv) {
	return Vector_V_sum_M_mul_V<MatrixN<T>, VectorN<T>>(v, Mv.M, Mv.a);
}

template <typename T> const Vector_V_sub_M_mul_V<MatrixN<T>, VectorN<T>> operator-(const VectorN<T> &v, const Vector_M_mul_V<MatrixN<T>, VectorN<T>> &Mv) {
	return Vector_V_sub_M_mul_V<MatrixN<T>, VectorN<T>>(v, Mv.M, Mv.a);
}

template <typename T> const Vector_V_sum_S_mul_M_mul_V<T, MatrixN<T>, VectorN<T>> operator+(const VectorN<T> &v, const Vector_S_mul_M_mul_V<T, MatrixN<T>, VectorN<T>> &sMv) {
	return Vector_V_sum_S_mul_M_mul_V<T, MatrixN<T>, VectorN<T>>(v, sMv.s, sMv.M, sMv.a);
}

template <typename T> const Vector_V_sub_S_mul_M_mul_V<T, MatrixN<T>, VectorN<T>> operator-(const VectorN<T> &v, const Vector_S_mul_M_mul_V<T, MatrixN<T>, VectorN<T>> &sMv) {
	return Vector_V_sub_S_mul_M_mul_V<T, MatrixN<T>, VectorN<T>>(v, sMv.s, sMv.M, sMv.a);
}

#undef __ASSERT

typedef MatrixN<float>			MatrixNf;
typedef MatrixN<double>			MatrixNd;

}; // namespace arithmetic

#endif