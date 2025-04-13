#ifndef _ARITHMETIC_MATRIX3_H_
#define _ARITHMETIC_MATRIX3_H_

#include "../../../Basis/Basis/BasicLibrary.h"
#include "../ArrayUtil/ArrayUtil3.h"
#include "../Vector/Vector3.h"

namespace arithmetic {

#ifdef _DEBUG
#define __ASSERT(c) assert(c)
#else
#define __ASSERT(c) {}
#endif

//////////////////////////////////////////////////////////////////////////
//	Matrix3 Class Definition.
//////////////////////////////////////////////////////////////////////////
template <typename T> class Matrix3 {
public:
	typedef unsigned int		index_type;
	typedef index_type*			index_pointer ;
	typedef const index_type*	index_const_pointer ;
	typedef index_type&			index_reference ;
	typedef const index_type&	index_const_reference ;

	typedef T					value_type;
	typedef value_type*			value_pointer ;       
	typedef const value_type*	value_const_pointer ; 
	typedef value_type&			value_reference ;     
	typedef const value_type&	value_const_reference;

	static const unsigned int getDimension() { return 3; }

	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	inline Matrix3();
	inline Matrix3( const Matrix3<T> &src );
	inline Matrix3( T e0, T e1, T e2, T e3,T e4, T e5, T e6, T e7, T e8 );
	inline Matrix3( const T *m);
	inline Matrix3( const T m[3][3] );
	inline Matrix3( const Vector3<T>& v0, const Vector3<T>& v1, const Vector3<T>& v2 );	

	//////////////////////////////////////////////////////////////////////////
	// Setter methods
	inline void set( T d );
	inline void set( index_type i, T d) ;
	inline void set( index_type i, index_type j, T d) ;
	inline void set( const Matrix3 &da );
	inline void set( T e0, T e1, T e2, T e3,T e4, T e5, T e6, T e7, T e8 );
	inline void set( const T *da );
	inline void set( const T da[3][3] );
	inline void set( const Vector3<T>& v0, const Vector3<T>& v1, const Vector3<T>& v2 );	
	inline Matrix3 &operator=( const Matrix3 &da );
	inline Matrix3 &operator=( const T *da );  
	inline Matrix3 &operator=( const T da[3][3] );  

	//////////////////////////////////////////////////////////////////////////
	// Getter methods	
	inline const T get( index_type i ) const 						{ return data[i];			};
	inline const T get( index_type i, index_type j ) const			{ return data[idx(i,j)];	};	
	inline const T operator()( index_type i ) const					{ return data[i];			};
	inline T&	   operator()( index_type i )						{ return data[i];			};
	inline const T operator()( index_type i, index_type j) const	{ return data[idx(i,j)];	};
	inline T&	   operator()( index_type i, index_type j)			{ return data[idx(i,j)];	};

	inline const Vector3<T> getRow( index_type i ) const;
	inline const Vector3<T> getColumn( index_type i)  const;

	inline const T det() const;
	inline const T determinant() const { return det(); }
	inline const T trace() const;

	inline Matrix3 &identity();
	inline Matrix3 &zero();

	inline const T* const getptr() const { return data; }

	//////////////////////////////////////////////////////////////////////////
	// Matrix Manipulation
	const bool isOrthoNormal() const;

	//////////////////////////////////////////////////////////////////////////
	// Matrix Operation
	inline Matrix3 &invert();
	const Matrix3 getInverse();
	inline Matrix3 &transpose();
	inline const Matrix3 getTranspose();
	inline Matrix3 &invertTranspose();
	inline const Matrix3 getInverseTranspose();

	//////////////////////////////////////////////////////////////////////////
	// Comparison operators
	inline const bool operator==(const Matrix3<T>& da) const;
	inline const bool operator!=(const Matrix3<T>& da) const;

	//////////////////////////////////////////////////////////////////////////
	// In place arithmetic
	inline Matrix3 &operator*=(T d);
	inline Matrix3 &operator/=(T d);

	inline Matrix3 &operator+=(const Matrix3 &da);
	inline Matrix3 &operator-=(const Matrix3 &da);
	inline Matrix3 &operator*=(const Matrix3 &da);

	//////////////////////////////////////////////////////////////////////////
	// Component Operators
	const index_type indexOfMaxDiagonalComponent() const;

	//////////////////////////////////////////////////////////////////////////
	// Friend methods
	template <typename T> friend const Matrix3<T> operator-( const Matrix3<T> &a );
	template <typename T> friend const Matrix3<T> operator+( const Matrix3<T> &a, const Matrix3<T> &b );
	template <typename T> friend const Matrix3<T> operator-( const Matrix3<T> &a, const Matrix3<T> &b );
	template <typename T> friend const Matrix3<T> operator*( const Matrix3<T> &a, const Matrix3<T> &b );

	template <typename T> friend const Vector3<T> operator*( const Matrix3<T> &a, const Vector3<T> &b );

	template <typename T> friend const Matrix3<T> operator+( const Matrix3<T> &a, const T b );
	template <typename T> friend const Matrix3<T> operator-( const Matrix3<T> &a, const T b );
	template <typename T> friend const Matrix3<T> operator*( const Matrix3<T> &a, const T b );
	template <typename T> friend const Matrix3<T> operator/( const Matrix3<T> &a, const T b );

	template <typename T> friend const Matrix3<T> operator+( const T a, const Matrix3<T> &b );
	template <typename T> friend const Matrix3<T> operator-( const T a, const Matrix3<T> &b );
	template <typename T> friend const Matrix3<T> operator*( const T a, const Matrix3<T> &b );

	template <typename T> friend std::ostream &operator<<( std::ostream &os, const Matrix3<T> &m );

private :
	inline const int idx(index_type row, index_type col) const { return (row + col*3); }

public : 
	T data[9];
};

//////////////////////////////////////////////////////////////////////////
// Constructor
template <typename T> inline Matrix3<T>::Matrix3() {
	identity();	
}
template <typename T> inline Matrix3<T>::Matrix3( const Matrix3<T> &src ) {
	set(src);	
}
template <typename T> inline Matrix3<T>::Matrix3( T e0, T e1, T e2, T e3,T e4, T e5, T e6, T e7,T e8 ) { 
	set(e0,e1,e2,e3,e4,e5,e6,e7,e8);
}
template <typename T> inline Matrix3<T>::Matrix3( const T *m ) { 
	set(m);
}
template <typename T> inline Matrix3<T>::Matrix3( const T m[3][3] ) { 
	set(m);
}
template <typename T> inline Matrix3<T>::Matrix3( const Vector3<T>& v0, const Vector3<T>& v1, const Vector3<T>& v2 ) { 
	set(v0,v1,v2);
}

//////////////////////////////////////////////////////////////////////////
// Setter methods
template <typename T> inline void Matrix3<T>::set( T d ) {
	for ( index_type i = 0; i < 9; ++i )
		data[i] = d;
}

template <typename T> inline void Matrix3<T>::set( index_type i, T d ) {
	__ASSERT( i < 9 );
	data[i] = d;
}

template <typename T> inline void Matrix3<T>::set( index_type i, index_type j, T d ) {
	__ASSERT( i < 3 && j < 3 );
	data[idx(i,j)] = d;
}

template <typename T> inline void Matrix3<T>::set( const Matrix3 &da ) {
	for ( index_type i = 0; i < 9; ++i ) 
		data[i] = da(i);
}

template <typename T> inline void Matrix3<T>::set( T e0, T e1, T e2, T e3,T e4, T e5, T e6, T e7, T e8 ) {
	data[0]=e0;		data[1]=e1;		data[2]=e2;		
	data[3]=e3;		data[4]=e4;		data[5]=e5;		
	data[6]=e6;		data[7]=e7;		data[8]=e8;		
}

template <typename T> inline void Matrix3<T>::set( const T *da ) {
	for ( index_type i = 0; i < 9; ++i )
		data[i] = da[i];
}

template <typename T> inline void Matrix3<T>::set( const T da[3][3] ) {
	for ( index_type i = 0; i < 3; ++i )
		for( index_type j = 0; j < 3; ++j)
			data[idx(i,j)] = da[i][j];
}

template <typename T> inline void Matrix3<T>::set( const Vector3<T>& v0, const Vector3<T>& v1, const Vector3<T>& v2 ) { 
	set(v0[0], v0[1], v0[2]
		v1[0], v1[1], v1[2]
		v2[0], v2[1], v2[2]);		
}

template <typename T> inline Matrix3<T> &Matrix3<T>::operator=( const Matrix3 &da ) {
	set(da);
	return (*this);
}
template <typename T> inline Matrix3<T> &Matrix3<T>::operator=( const T *da ) {
	set(da);
	return (*this);
}
template <typename T> inline Matrix3<T> &Matrix3<T>::operator=( const T da[3][3] ) {
	set(da);
	return (*this);
}

//////////////////////////////////////////////////////////////////////////
// Getter methods
template <typename T> inline const Vector3<T> Matrix3<T>::getRow(index_type i) const {
	__ASSERT( i < 3 );
	if(i==0) return Vector3<T>(data[0], data[3], data[6]);
	if(i==1) return Vector3<T>(data[1], data[4], data[7]);
	if(i==2) return Vector3<T>(data[2], data[5], data[8]);
	assert(false);
	return Vector3<T>(data[0], data[3], data[6]);
}

template <typename T> inline const Vector3<T> Matrix3<T>::getColumn(index_type i)	const {
	__ASSERT( i < 3 );
	if(i==0) return Vector3<T>(data[0], data[1], data[2]);
	if(i==1) return Vector3<T>(data[3], data[4], data[5]);
	if(i==2) return Vector3<T>(data[6], data[7], data[8]);
	assert(false);
	return Vector3<T>(data[0], data[1], data[2]);	
}

template <typename T> inline const T Matrix3<T>::det() const {
	return data[0]*(data[4]*data[8]-data[7]*data[5])
		  -data[3]*(data[1]*data[8]-data[7]*data[2]) 
		  +data[6]*(data[1]*data[5]-data[4]*data[2]);
}

template <typename T> inline const T Matrix3<T>::trace() const {
	return (*this)(0,0) + (*this)(1,1) + (*this)(2,2);
}

template <typename T> inline Matrix3<T> &Matrix3<T>::identity()	{
	set(T(1.0), T(0.0), T(0.0), T(0.0), T(1.0), T(0.0), T(0.0), T(0.0), T(1.0));
	return (*this);
}

template <typename T> inline Matrix3<T> &Matrix3<T>::zero()	{
	set(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0));
	return (*this);		
}

//////////////////////////////////////////////////////////////////////////
// Matrix Operation
template <typename T> inline const bool Matrix3<T>::isOrthoNormal() const {
	Vector3<T> x = getRow(0), y =  getRow(1), z = getRow(2);
	if (!IsAlmostSame(x.sqrmag(),(T)1)) return false;
	if (!IsAlmostSame(y.sqrmag(),(T)1)) return false;
	if (!IsAlmostSame(z.sqrmag(),(T)1)) return false;
	if (!IsAlmostZero(dot(x,y))) return false;
	if (!IsAlmostZero(dot(y,z))) return false;
	if (!IsAlmostZero(dot(z,x))) return false;
	
	return true;
}

//////////////////////////////////////////////////////////////////////////
// Matrix Operation
template <typename T> inline Matrix3<T> &Matrix3<T>::invert()	{
	return (*this)=getInverse();
}

template <typename T> const Matrix3<T> Matrix3<T>::getInverse()	{

	T invDet = 1.f/det();
	return Matrix3<T>(	(data[4]*data[8]-data[5]*data[7])*invDet,
						(data[2]*data[7]-data[1]*data[8])*invDet,
						(data[1]*data[5]-data[2]*data[4])*invDet,
						(data[5]*data[6]-data[3]*data[8])*invDet,
						(data[0]*data[8]-data[2]*data[6])*invDet,
						(data[2]*data[3]-data[0]*data[5])*invDet,
						(data[3]*data[7]-data[4]*data[6])*invDet,
						(data[1]*data[6]-data[0]*data[7])*invDet,
						(data[0]*data[4]-data[1]*data[3])*invDet );
		
}

template <typename T> inline Matrix3<T> &Matrix3<T>::transpose()	{
	return (*this)=getTranspose();
}

template <typename T> inline const Matrix3<T> Matrix3<T>::getTranspose()	{
	return Matrix3<T>(data[0], data[3], data[6], 
					  data[1], data[4], data[7], 
					  data[2], data[5], data[8] );
}

template <typename T> inline Matrix3<T> &Matrix3<T>::invertTranspose()	{
	return (*this)=getInverseTranspose();
}

template <typename T> inline const Matrix3<T> Matrix3<T>::getInverseTranspose()	{
	Matrix3<T> result = getInverse();
	return result.getTranspose();
}

//////////////////////////////////////////////////////////////////////////
// Comparison operators
template <typename T> inline const bool Matrix3<T>::operator==(const Matrix3<T>& da) const	{

	for(int i=0; i<9; ++i)	{
		if(data[i] != da.data[i])
			return false;
	}
	return true;
}

template <typename T> inline const bool Matrix3<T>::operator!=(const Matrix3<T>& da) const	{
	return !((*this)==da);
}

//////////////////////////////////////////////////////////////////////////
// In place arithmetic
template <typename T> inline Matrix3<T> &Matrix3<T>::operator*=(T d)	{
	return ((*this) = (*this) * d);
}

template <typename T> inline Matrix3<T> &Matrix3<T>::operator/=(T d)	{
	return ((*this) = (*this) / d);
}

template <typename T> inline Matrix3<T> &Matrix3<T>::operator+=(const Matrix3 &da)	{
	return ((*this) = (*this) + da);
}

template <typename T> inline Matrix3<T> &Matrix3<T>::operator-=(const Matrix3 &da)	{
	return ((*this) = (*this) - da);
}

template <typename T> inline Matrix3<T> &Matrix3<T>::operator*=(const Matrix3 &da)	{
	return ((*this) = (*this) * da);
}

//////////////////////////////////////////////////////////////////////////
// Component Operator
template <typename T> 
const typename Matrix3<T>::index_type Matrix3<T>::indexOfMaxDiagonalComponent() const {
	Vector3<T> diag((*this)(0,0), (*this)(1,1), (*this)(2,2));
	return diag.indexOfMaxComponent();
}

//////////////////////////////////////////////////////////////////////////
// Friend methods
template <typename T> inline const Matrix3<T> operator-( const Matrix3<T> &a ) {
	return Matrix3<T>(-a.data[0],-a.data[1],-a.data[2],
					  -a.data[3],-a.data[4],-a.data[5],
					  -a.data[6],-a.data[7],-a.data[8] );
}
template <typename T> inline const Matrix3<T> operator+( const Matrix3<T> &a, const Matrix3<T> &b ) {
	return Matrix3<T>(a.data[0]+b.data[0],a.data[1]+b.data[1],a.data[2]+b.data[2],
					  a.data[3]+b.data[3],a.data[4]+b.data[4],a.data[5]+b.data[5],
					  a.data[6]+b.data[6],a.data[7]+b.data[7],a.data[8]+b.data[8] );
}
template <typename T> inline const Matrix3<T> operator-( const Matrix3<T> &a, const Matrix3<T> &b ) {
	return Matrix3<T>(a.data[0]-b.data[0],a.data[1]-b.data[1],a.data[2]-b.data[2],
					  a.data[3]-b.data[3],a.data[4]-b.data[4],a.data[5]-b.data[5],
					  a.data[6]-b.data[6],a.data[7]-b.data[7],a.data[8]-b.data[8] );
}
template <typename T> inline const Matrix3<T> operator*( const Matrix3<T> &a, const Matrix3<T> &b ) {
	return Matrix3<T>( a.data[0]*b.data[0] + a.data[3]*b.data[1] + a.data[6]*b.data[2],
					   a.data[1]*b.data[0] + a.data[4]*b.data[1] + a.data[7]*b.data[2],
					   a.data[2]*b.data[0] + a.data[5]*b.data[1] + a.data[8]*b.data[2],					   
					   
					   a.data[0]*b.data[3] + a.data[3]*b.data[4] + a.data[6]*b.data[5],					   
					   a.data[1]*b.data[3] + a.data[4]*b.data[4] + a.data[7]*b.data[5],
					   a.data[2]*b.data[3] + a.data[5]*b.data[4] + a.data[8]*b.data[5],					   
					   
					   a.data[0]*b.data[6] + a.data[3]*b.data[7] + a.data[6]*b.data[8],					   
					   a.data[1]*b.data[6] + a.data[4]*b.data[7] + a.data[7]*b.data[8],					   
					   a.data[2]*b.data[6] + a.data[5]*b.data[7] + a.data[8]*b.data[8] );
}

template <typename T> inline const Vector3<T> operator*( const Matrix3<T> &a, const Vector3<T> &b ) {
	return Vector3<T>( a.data[0]*b[0] + a.data[3]*b[1] + a.data[6]*b[2],
					   a.data[1]*b[0] + a.data[4]*b[1] + a.data[7]*b[2],
					   a.data[2]*b[0] + a.data[5]*b[1] + a.data[8]*b[2] );
}

template <typename T> inline const Matrix3<T> operator+( const Matrix3<T> &a, const T b ) {
	return Matrix3<T>( a.data[0]+b,	a.data[1]+b, a.data[2]+b, 
					   a.data[3]+b, a.data[4]+b, a.data[5]+b, 
					   a.data[6]+b, a.data[7]+b, a.data[8]+b	);
}
template <typename T> inline const Matrix3<T> operator-( const Matrix3<T> &a, const T b ) {
	return Matrix3<T>( a.data[0]-b,	a.data[1]-b, a.data[2]-b, 
					   a.data[3]-b, a.data[4]-b, a.data[5]-b, 
					   a.data[6]-b, a.data[7]-b, a.data[8]-b	);
}
template <typename T> inline const Matrix3<T> operator*( const Matrix3<T> &a, const T b ) {
	return Matrix3<T>( a.data[0]*b,	a.data[1]*b, a.data[2]*b, 
					   a.data[3]*b, a.data[4]*b, a.data[5]*b, 
					   a.data[6]*b, a.data[7]*b, a.data[8]*b	);
}
template <typename T> inline const Matrix3<T> operator/( const Matrix3<T> &a, const T b ) {
	return Matrix3<T>( a.data[0]/b,	a.data[1]/b, a.data[2]/b, 
					   a.data[3]/b, a.data[4]/b, a.data[5]/b, 
					   a.data[6]/b, a.data[7]/b, a.data[8]/b	);
}

template <typename T> inline const Matrix3<T> operator+( T a, const Matrix3<T> &b ) {
	return (b+a);
}
template <typename T> inline const Matrix3<T> operator-( T a, const Matrix3<T> &b ) {
	return (-b+a);
}
template <typename T> inline const Matrix3<T> operator*( T a, const Matrix3<T> &b ) {
	return (b*a);
}

template <typename T> std::ostream &operator<<( std::ostream &os, const Matrix3<T> &m ) {
	os	<< "[" << m.data[0] << ", " << m.data[3] << ", " << m.data[6] << std::endl
		<< " " << m.data[1] << ", " << m.data[4] << ", " << m.data[7] << std::endl
		<< " " << m.data[2] << ", " << m.data[5] << ", " << m.data[8] << "]";
	return os;
}

template <typename T>
inline Matrix3<T> Vector3<T>::getDualMatrix() const {
	return Matrix3<T>(0.0f, data[2], -data[1],	
					 -data[2], 0.0f, data[0],
					  data[1], -data[0], 0.0f);			// colume-wise order
}

template <typename T>
inline Matrix3<T> Vector3<T>::getTensorProduct() const {
	Matrix3<T> mat;
	getTensorProduct(mat);
	return mat;
}

template <typename T>
inline void Vector3<T>::getTensorProduct(Matrix3<T>& mat) const {
	for(index_type i=0;i<3;++i) 
		for(index_type j=0;j<3;++j) 
			mat(i,j) = data[i] * data[j];
}

typedef Matrix3<float>			Matrix3f;
typedef Matrix3<double>			Matrix3d;

#undef __ASSERT

};		// namespace arithmetic

#endif