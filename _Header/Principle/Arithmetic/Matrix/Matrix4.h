#ifndef _ARITHMETIC_MATRIX4_H_
#define _ARITHMETIC_MATRIX4_H_

#include "../../../Basis/Basis/BasicLibrary.h"
#include "../ArrayUtil/ArrayUtil3.h"
#include "../Vector/Vector3.h"
#include "../Vector/Vector4.h"

namespace arithmetic {

using namespace basis;

#ifdef _DEBUG
#define __ASSERT(c) assert(c)
#else
#define __ASSERT(c) {}
#endif

//////////////////////////////////////////////////////////////////////////
//	Matrix4 Class Definition.
//////////////////////////////////////////////////////////////////////////
template <typename T> class Matrix4 {
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

	static const unsigned int getDimension() { return 4; }

	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	inline Matrix4();
	template<typename S> inline Matrix4( const Matrix4<S> &src );
	template<typename S> explicit inline Matrix4( const Matrix3<S> &src );
	inline Matrix4( T e0, T e1, T e2, T e3,T e4, T e5, T e6, T e7, T e8, T e9, T e10, T e11,T e12, T e13, T e14, T e15);
	inline Matrix4( const T *m);
	template<typename S> inline Matrix4( const S m[4][4] );
	inline Matrix4( const Vector4<T>& v0, const Vector4<T>& v1, const Vector4<T>& v2, const Vector4<T>& v3 );	

	//////////////////////////////////////////////////////////////////////////
	// Setter methods
	inline void set( T d );
	inline void set( index_type i, T d) ;
	inline void set( index_type i, index_type j, T d) ;
	template<typename S> inline void set( const Matrix4<S> &da );
	template<typename S> inline void set( const Matrix3<S> &da );
	inline void set( T e0, T e1, T e2, T e3,T e4, T e5, T e6, T e7,T e8, T e9, T e10, T e11,T e12, T e13, T e14, T e15);
	inline void set( const T *da );
	template<typename S> inline void set( const S da[4][4] );
	inline void set( const Vector4<T>& v0, const Vector4<T>& v1, const Vector4<T>& v2, const Vector4<T>& v3 );
	template<typename S> inline Matrix4 &operator=( const Matrix4<S> &da );
	template<typename S> inline Matrix4 &operator=( const Matrix3<S> &da );
	inline Matrix4 &operator=( const T *da );  
	template<typename S> inline Matrix4 &operator=( const S da[4][4] );  

	//////////////////////////////////////////////////////////////////////////
	// Getter methods
	inline const T get( index_type i ) const 						{ return data[i];			};
	inline const T get( index_type i, index_type j ) const			{ return data[idx(i,j)];	};	
	inline const T operator()( index_type i ) const					{ return data[i];			};
	inline T&	   operator()( index_type i )						{ return data[i];			};
	inline const T operator()( index_type i, index_type j) const	{ return data[idx(i,j)];	};
	inline T&	   operator()( index_type i, index_type j)			{ return data[idx(i,j)];	};

	inline const Vector4<T> getRow( index_type i ) const;
	inline const Vector4<T> getColumn( index_type i)  const;

	inline Matrix4 &identity();
	inline Matrix4 &zero();

	inline const T* const getptr() const { return data; }

	//////////////////////////////////////////////////////////////////////////
	// Matrix Operation
	inline Matrix4 &invert();
	inline const Matrix4 getInverse();
	inline Matrix4 &transpose();
	inline const Matrix4 getTranspose();
	inline Matrix4 &invertTranspose();
	const Matrix4 getInverseTranspose();

	inline Matrix4 &affineInvert();
	inline const Matrix4 getAffineInverse() const;
	inline Matrix4 &affineInvertTranspose();
	inline const Matrix4 getAffineInverseTranspose() const;

	//////////////////////////////////////////////////////////////////////////
	// Vector Operation
	// Rotate a 3d vector by rotation part
	void rotateVector3(Vector3<T> &v) const			{ v = getRotatedVector3(v);			}
	void inverseRotateVector3(Vector3<T> & v) const	{ v = getInverseRotatedVector3(v);	}

	inline const Vector3<T> getRotatedVector3(const Vector3<T> & v) const;
	inline const Vector3<T> getInverseRotatedVector3(const Vector3<T> & v) const;

	// Translate a 3d vector by translation part
	void translateVector3(Vector3<T> & v) const			{ v = getTranslatedVector3(v);			}
	void inverseTranslateVector3(Vector3<T> & v) const	{ v = getInverseTranslatedVector3(v);	}

	inline const Vector3<T> getTranslatedVector3(const Vector3<T> & v) const;
	inline const Vector3<T> getInverseTranslatedVector3(const Vector3<T> & v) const;

	//////////////////////////////////////////////////////////////////////////
	// Special Matrix Setting
	inline const Matrix4 &setTranslation(const Vector3<T>& translation);
	inline const Matrix4 &setScale(const Vector3<T>& scaleFactor);
	inline const Matrix4 &setUniformScale(const T scaleFactor);
	const Matrix4 &setRotationAxis(const T angle, const Vector3<T>& axis);
	const Matrix4 &setRotationX(const T angle);
	const Matrix4 &setRotationY(const T angle);
	const Matrix4 &setRotationZ(const T angle);
	const Matrix4 &setRotationEuler(const T angleX, const T angleY, const T angleZ);

	// Set parts of the matrix
	inline Matrix4 &setTranslationPart(const Vector3<T> & translation);
	Matrix4 &setRotationPartEuler(const T angleX, const T angleY, const T angleZ);
	Matrix4 &setRotationPartEuler(const Vector3<T> & rotations)	{
		return setRotationPartEuler(rotations.x, rotations.y, rotations.z);
	}

	// Set View Matrix
	const Matrix4 &setPerspective(T left, T right, T bottom, T top, T n, T f);
	const Matrix4 &setPerspective(T fovy, T aspect, T n, T f);
	const Matrix4 &setOrtho(T left, T right, T bottom, T top, T n, T f);
	const Matrix4 &setView(const Vector3<T> &eye, const Vector3<T> &center, const Vector3<T> &up);
	const Matrix4 &setView(const T eye[3], const T center[3], const T up[3]);

	//////////////////////////////////////////////////////////////////////////
	// Comparison operators
	inline const bool operator==(const Matrix4<T>& da) const;
	inline const bool operator!=(const Matrix4<T>& da) const;

	//////////////////////////////////////////////////////////////////////////
	// In place arithmetic
	inline Matrix4 &operator*=(T d);
	inline Matrix4 &operator/=(T d);

	inline Matrix4 &operator+=(const Matrix4 &da);
	inline Matrix4 &operator-=(const Matrix4 &da);
	inline Matrix4 &operator*=(const Matrix4 &da);

	//////////////////////////////////////////////////////////////////////////
	// Friend methods
	template <typename T> friend const Matrix4<T> operator-( const Matrix4<T> &a );
	template <typename T> friend const Matrix4<T> operator+( const Matrix4<T> &a, const Matrix4<T> &b );
	template <typename T> friend const Matrix4<T> operator-( const Matrix4<T> &a, const Matrix4<T> &b );
	template <typename T> friend const Matrix4<T> operator*( const Matrix4<T> &a, const Matrix4<T> &b );

	template <typename T> friend const Vector4<T> operator*( const Matrix4<T> &a, const Vector4<T> &b );
	template <typename T> friend const Vector3<T> operator*( const Matrix4<T> &a, const Vector3<T> &b );

	template <typename T> friend const Matrix4<T> operator+( const Matrix4<T> &a, const T b );
	template <typename T> friend const Matrix4<T> operator-( const Matrix4<T> &a, const T b );
	template <typename T> friend const Matrix4<T> operator*( const Matrix4<T> &a, const T b );
	template <typename T> friend const Matrix4<T> operator/( const Matrix4<T> &a, const T b );

	template <typename T> friend const Matrix4<T> operator+( const T a, const Matrix4<T> &b );
	template <typename T> friend const Matrix4<T> operator-( const T a, const Matrix4<T> &b );
	template <typename T> friend const Matrix4<T> operator*( const T a, const Matrix4<T> &b );

	template <typename T> friend std::ostream &operator<<( std::ostream &strm, const Matrix4<T> &m );

private :
	inline const index_type idx(index_type row, index_type col) const { return (row + col*4); }

public : 
	T data[16];
};

//////////////////////////////////////////////////////////////////////////
// Constructor
template <typename T> inline Matrix4<T>::Matrix4() {
	identity();	
}
template <typename T> template<typename S> inline Matrix4<T>::Matrix4( const Matrix4<S> &src ) {
	set(src);	
}
template <typename T> template<typename S> inline Matrix4<T>::Matrix4( const Matrix3<S> &src ) {
	set(src);	
}
template <typename T> inline Matrix4<T>::Matrix4( T e0, T e1, T e2, T e3,T e4, T e5, T e6, T e7,T e8, T e9, T e10, T e11,T e12, T e13, T e14, T e15) { 
	set(e0,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12,e13,e14,e15);
}
template <typename T> inline Matrix4<T>::Matrix4( const T *m ) { 
	set(m);
}
template <typename T> template<typename S> inline Matrix4<T>::Matrix4( const S m[4][4] ) { 
	set(m);
}
template <typename T> inline Matrix4<T>::Matrix4( const Vector4<T>& v0, const Vector4<T>& v1, const Vector4<T>& v2, const Vector4<T>& v3 ) { 
	set(v0,v1,v2,v3);
}

//////////////////////////////////////////////////////////////////////////
// Setter methods
template <typename T> inline void Matrix4<T>::set( T d ) {
	for ( index_type i = 0; i < 16; ++i )
		data[i] = d;
}

template <typename T> inline void Matrix4<T>::set( index_type i, T d ) {
	__ASSERT( i < 16 );
	data[i] = d;
}

template <typename T> inline void Matrix4<T>::set( index_type i, index_type j, T d ) {
	__ASSERT( i < 4 && j < 4 );
	data[idx(i,j)] = d;
}

template <typename T> template<typename S> inline void Matrix4<T>::set( const Matrix4<S> &da ) {
	for ( index_type i = 0; i < 16; ++i ) 
		data[i] = da(i);
}

template <typename T> template<typename S> inline void Matrix4<T>::set( const Matrix3<S> &da ) {
	identity();
	for ( index_type i = 0; i < 3; ++i ) {
		for( index_type j = 0; j < 3; ++j ) {
			(*this)(i,j) = da(i,j);
		}
	}		
}

template <typename T> inline void Matrix4<T>::set( T e0, T e1, T e2, T e3,T e4, T e5, T e6, T e7,T e8, T e9, T e10, T e11,T e12, T e13, T e14, T e15) {
	data[0]=e0;		data[1]=e1;		data[2]=e2;		data[3]=e3;
	data[4]=e4;		data[5]=e5;		data[6]=e6;		data[7]=e7;
	data[8]=e8;		data[9]=e9;		data[10]=e10;	data[11]=e11;
	data[12]=e12;	data[13]=e13;	data[14]=e14;	data[15]=e15;
}

template <typename T> inline void Matrix4<T>::set( const T *da ) {
	for ( index_type i = 0; i < 16; ++i )
		data[i] = da[i];
}

template <typename T> template<typename S> inline void Matrix4<T>::set( const S da[4][4] ) {
	for ( index_type i = 0; i < 4; ++i )
		for( index_type j = 0; j < 4; ++j)
			data[idx(i,j)] = da[i][j];
}

template <typename T> inline void Matrix4<T>::set( const Vector4<T>& v0, const Vector4<T>& v1, const Vector4<T>& v2, const Vector4<T> &v3 ) { 
	set(v0[0], v0[1], v0[2], v0[3], 
		v1[0], v1[1], v1[2], v1[3], 
		v2[0], v2[1], v2[2], v2[3], 
		v3[0], v3[1], v3[2], v3[3]);
}

template <typename T> template<typename S> inline Matrix4<T> &Matrix4<T>::operator=( const Matrix4<S> &da ) {
	set(da);
	return (*this);
}
template <typename T> template<typename S> inline Matrix4<T> &Matrix4<T>::operator=( const Matrix3<S> &da ) {
	set(da);
	return (*this);
}
template <typename T> inline Matrix4<T> &Matrix4<T>::operator=( const T *da ) {
	set(da);
	return (*this);
}
template <typename T> template<typename S> inline Matrix4<T> &Matrix4<T>::operator=( const S da[4][4] ) {
	set(da);
	return (*this);
}

//////////////////////////////////////////////////////////////////////////
// Getter methods
template <typename T> inline const Vector4<T> Matrix4<T>::getRow(index_type i) const {
	__ASSERT( i < 4 );
	if(i==0) return Vector4<T>(data[0], data[4], data[8], data[12]);
	if(i==1) return Vector4<T>(data[1], data[5], data[9], data[13]);
	if(i==2) return Vector4<T>(data[2], data[6], data[10], data[14]);
	if(i==3) return Vector4<T>(data[3], data[7], data[11], data[15]);
}

template <typename T> inline const Vector4<T> Matrix4<T>::getColumn(index_type i)	const {
	__ASSERT( i < 4 );
	if(i==0) return Vector4<T>(data[0], data[1], data[2], data[3]);
	if(i==1) return Vector4<T>(data[4], data[5], data[6], data[7]);
	if(i==2) return Vector4<T>(data[8], data[9], data[10], data[11]);
	if(i==3) return Vector4<T>(data[12], data[13], data[14], data[15]);
}

template <typename T> inline Matrix4<T> &Matrix4<T>::identity()	{
	set(T(1.0), T(0.0), T(0.0), T(0.0), T(0.0), T(1.0), T(0.0), T(0.0), T(0.0), T(0.0), T(1.0), T(0.0), T(0.0), T(0.0), T(0.0), T(1.0));
	return (*this);
}

template <typename T> inline Matrix4<T> &Matrix4<T>::zero()	{
	set(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0));
	return (*this);
}

//////////////////////////////////////////////////////////////////////////
// Matrix Operation
template <typename T> inline Matrix4<T> &Matrix4<T>::invert()	{
	return (*this)=getInverse();
}

template <typename T> inline const Matrix4<T> Matrix4<T>::getInverse()	{
	Matrix4<T> result = getInverseTranspose();
	return result.getTranspose();
}

template <typename T> inline Matrix4<T> &Matrix4<T>::transpose()	{
	return (*this)=getTranspose();
}

template <typename T> inline const Matrix4<T> Matrix4<T>::getTranspose()	{
	return Matrix4<T>(	data[ 0], data[ 4], data[ 8], data[12],
						data[ 1], data[ 5], data[ 9], data[13],
						data[ 2], data[ 6], data[10], data[14],
						data[ 3], data[ 7], data[11], data[15]);
}

template <typename T> inline Matrix4<T> &Matrix4<T>::invertTranspose()	{
	return (*this)=getInverseTranspose();
}

template <typename T> const Matrix4<T> Matrix4<T>::getInverseTranspose()	{
	
	Matrix4<T> result;

	T tmp[12];												//temporary pair storage
	T det;													//determinant

	//calculate pairs for first 8 elements (cofactors)
	tmp[0] = data[10] * data[15];
	tmp[1] = data[11] * data[14];
	tmp[2] = data[9] * data[15];
	tmp[3] = data[11] * data[13];
	tmp[4] = data[9] * data[14];
	tmp[5] = data[10] * data[13];
	tmp[6] = data[8] * data[15];
	tmp[7] = data[11] * data[12];
	tmp[8] = data[8] * data[14];
	tmp[9] = data[10] * data[12];
	tmp[10] = data[8] * data[13];
	tmp[11] = data[9] * data[12];

	//calculate first 8 elements (cofactors)
	result.set(0, tmp[0]*data[5] + tmp[3]*data[6] + tmp[4]*data[7] -	tmp[1]*data[5] - tmp[2]*data[6] - tmp[5]*data[7]);
	result.set(1, tmp[1]*data[4] + tmp[6]*data[6] + tmp[9]*data[7] -	tmp[0]*data[4] - tmp[7]*data[6] - tmp[8]*data[7]);
	result.set(2, tmp[2]*data[4] + tmp[7]*data[5] + tmp[10]*data[7] -	tmp[3]*data[4] - tmp[6]*data[5] - tmp[11]*data[7]);
	result.set(3, tmp[5]*data[4] + tmp[8]*data[5] + tmp[11]*data[6] -	tmp[4]*data[4] - tmp[9]*data[5] - tmp[10]*data[6]);
	result.set(4, tmp[1]*data[1] + tmp[2]*data[2] + tmp[5]*data[3] -	tmp[0]*data[1] - tmp[3]*data[2] - tmp[4]*data[3]);
	result.set(5, tmp[0]*data[0] + tmp[7]*data[2] + tmp[8]*data[3] -	tmp[1]*data[0] - tmp[6]*data[2] - tmp[9]*data[3]);
	result.set(6, tmp[3]*data[0] + tmp[6]*data[1] + tmp[11]*data[3] -	tmp[2]*data[0] - tmp[7]*data[1] - tmp[10]*data[3]);
	result.set(7, tmp[4]*data[0] + tmp[9]*data[1] + tmp[10]*data[2] -	tmp[5]*data[0] - tmp[8]*data[1] - tmp[11]*data[2]);

	//calculate pairs for second 8 elements (cofactors)
	tmp[0] = data[2]*data[7];
	tmp[1] = data[3]*data[6];
	tmp[2] = data[1]*data[7];
	tmp[3] = data[3]*data[5];
	tmp[4] = data[1]*data[6];
	tmp[5] = data[2]*data[5];
	tmp[6] = data[0]*data[7];
	tmp[7] = data[3]*data[4];
	tmp[8] = data[0]*data[6];
	tmp[9] = data[2]*data[4];
	tmp[10] = data[0]*data[5];
	tmp[11] = data[1]*data[4];

	//calculate second 8 elements (cofactors)
	result.set(8,		tmp[0]*data[13] + tmp[3]*data[14] + tmp[4]*data[15]	-	tmp[1]*data[13] - tmp[2]*data[14] - tmp[5]*data[15]);
	result.set(9,		tmp[1]*data[12] + tmp[6]*data[14] + tmp[9]*data[15]	-	tmp[0]*data[12] - tmp[7]*data[14] - tmp[8]*data[15]);
	result.set(10,		tmp[2]*data[12] + tmp[7]*data[13] + tmp[10]*data[15] -	tmp[3]*data[12] - tmp[6]*data[13] - tmp[11]*data[15]);
	result.set(11,		tmp[5]*data[12] + tmp[8]*data[13] + tmp[11]*data[14] -	tmp[4]*data[12] - tmp[9]*data[13] - tmp[10]*data[14]);
	result.set(12,		tmp[2]*data[10] + tmp[5]*data[11] + tmp[1]*data[9]	-	tmp[4]*data[11] - tmp[0]*data[9] - tmp[3]*data[10]);
	result.set(13,		tmp[8]*data[11] + tmp[0]*data[8] + tmp[7]*data[10]	-	tmp[6]*data[10] - tmp[9]*data[11] - tmp[1]*data[8]);
	result.set(14,		tmp[6]*data[9] + tmp[11]*data[11] + tmp[3]*data[8]	-	tmp[10]*data[11] - tmp[2]*data[8] - tmp[7]*data[9]);
	result.set(15,		tmp[10]*data[10] + tmp[4]*data[8] + tmp[9]*data[9]	-	tmp[8]*data[9] - tmp[11]*data[10] - tmp[5]*data[8]);

	// calculate determinant
	det	= data[0]*result.get(0) + data[1]*result.get(1) + data[2]*result.get(2) + data[3]*result.get(3);

	if(IsAlmostZero(det))	{
		Matrix4<T> id;
		return id;
	}

	return result / det;
}

//Invert if only composed of rotations & translations
template <typename T> inline Matrix4<T>& Matrix4<T>::affineInvert()	{
	return (*this) = getAffineInverse();
}

template <typename T> inline const Matrix4<T> Matrix4<T>::getAffineInverse() const	{
	//return the transpose of the rotation part
	//and the negative of the inverse rotated translation part
	return Matrix4<T>(	data[0],	data[4],	data[8],	0.0f,
						data[1],	data[5],	data[9],	0.0f,
						data[2],	data[6],	data[10],	0.0f,
						-(data[0]*data[12]+data[1]*data[13]+data[2]*data[14]),
						-(data[4]*data[12]+data[5]*data[13]+data[6]*data[14]),
						-(data[8]*data[12]+data[9]*data[13]+data[10]*data[14]),
						1.0f);
}

template <typename T> inline Matrix4<T>& Matrix4<T>::affineInvertTranspose()	{
	return (*this)=getAffineInverseTranspose();
}

template <typename T> inline const Matrix4<T> Matrix4<T>::getAffineInverseTranspose() const	{
	//return the transpose of the rotation part
	//and the negative of the inverse rotated translation part
	//transposed
	return Matrix4<T>(	data[0],	data[1],	data[2],	-(data[0]*data[12]+data[1]*data[13]+data[2]*data[14]),
						data[4],	data[5],	data[6],	-(data[4]*data[12]+data[5]*data[13]+data[6]*data[14]),
						data[8],	data[9],	data[10],	-(data[8]*data[12]+data[9]*data[13]+data[10]*data[14]),
						0.0f, 0.0f, 0.0f, 1.0f);
}

//////////////////////////////////////////////////////////////////////////
// Vector Operation
// Rotate a 3d vector by rotation part
template <typename T> inline const Vector3<T> Matrix4<T>::getRotatedVector3(const Vector3<T> & v) const	{
	return Vector3<T>(data[0]*v[0] + data[4]*v[1] + data[8]*v[2],
					  data[1]*v[0] + data[5]*v[1] + data[9]*v[2],
					  data[2]*v[0] + data[6]*v[1] + data[10]*v[2]);
}

template <typename T> inline const Vector3<T> Matrix4<T>::getInverseRotatedVector3(const Vector3<T> & v) const	{
	return Vector3<T>(data[0]*v[0] + data[1]*v[1] + data[2]*v[2],
					  data[4]*v[0] + data[5]*v[1] + data[6]*v[2],
					  data[8]*v[0] + data[9]*v[1] + data[10]*v[2]);
}

// Translate a 3d vector by translation part
template <typename T> inline const Vector3<T> Matrix4<T>::getTranslatedVector3(const Vector3<T> & v) const	{
	return Vector3<T>(v[0]+data[12], v[1]+data[13], v[2]+data[14]);
}

template <typename T> inline const Vector3<T> Matrix4<T>::getInverseTranslatedVector3(const Vector3<T> & v) const	{
	return Vector3<T>(v[0]-data[12], v[1]-data[13], v[2]-data[14]);
}

//////////////////////////////////////////////////////////////////////////
// Special Matrix Setting
template <typename T> inline const Matrix4<T> &Matrix4<T>::setTranslation(const Vector3<T>& translation)	{
	identity();
	setTranslationPart(translation);	
	return (*this);
}

template <typename T> inline const Matrix4<T> &Matrix4<T>::setScale(const Vector3<T>& scaleFactor)	{
	identity();
	data[0] = scaleFactor[0];
	data[5] = scaleFactor[1];
	data[10] = scaleFactor[2];
	return (*this);
}

template <typename T> inline const Matrix4<T> &Matrix4<T>::setUniformScale(const T scaleFactor)	{
	identity();
	data[0]=data[5]=data[10]=scaleFactor;
	return (*this);
}

template <typename T> const Matrix4<T> &Matrix4<T>::setRotationAxis(const T angle, const Vector3<T>& axis)	{
	Vector3<T> u = unitVector(axis);

	T sinAngle = sin(static_cast<T>(basis::PI) * angle / 180.0f);
	T cosAngle = cos(static_cast<T>(basis::PI) * angle / 180.0f);
	T oneMinusCosAngle = 1.0f - cosAngle;

	identity();

	data[0]=(u[0])*(u[0]) + cosAngle*(1-(u[0])*(u[0]));
	data[4]=(u[0])*(u[1])*(oneMinusCosAngle) - sinAngle*u[2];
	data[8]=(u[0])*(u[2])*(oneMinusCosAngle) + sinAngle*u[1];

	data[1]=(u[0])*(u[1])*(oneMinusCosAngle) + sinAngle*u[2];
	data[5]=(u[1])*(u[1]) + cosAngle*(1-(u[1])*(u[1]));
	data[9]=(u[1])*(u[2])*(oneMinusCosAngle) - sinAngle*u[0];

	data[2]=(u[0])*(u[2])*(oneMinusCosAngle) - sinAngle*u[1];
	data[6]=(u[1])*(u[2])*(oneMinusCosAngle) + sinAngle*u[0];
	data[10]=(u[2])*(u[2]) + cosAngle*(1-(u[2])*(u[2]));

	return (*this);
}

template <typename T> const Matrix4<T> &Matrix4<T>::setRotationX(const T angle)	{
	identity();

	data[5]=cos(static_cast<T>(basis::PI)*angle/180);
	data[6]=sin(static_cast<T>(basis::PI)*angle/180);

	data[9]=-data[6];
	data[10]=data[5];

	return (*this);
}

template <typename T> const Matrix4<T> &Matrix4<T>::setRotationY(const T angle)	{
	identity();

	data[0]=cos(static_cast<T>(basis::PI)*angle/180);
	data[2]=-sin(static_cast<T>(basis::PI)*angle/180);

	data[8]=-data[2];
	data[10]=data[0];

	return (*this);
}

template <typename T> const Matrix4<T> &Matrix4<T>::setRotationZ(const T angle)	{
	identity();

	data[0]=cos(static_cast<T>(basis::PI)*angle/180);
	data[1]=sin(static_cast<T>(basis::PI)*angle/180);

	data[4]=-data[1];
	data[5]=data[0];

	return (*this);
}

template <typename T> const Matrix4<T> &Matrix4<T>::setRotationEuler(const T angleX, const T angleY, const T angleZ)	{
	identity();
	return setRotationPartEuler(angleX, angleY, angleZ);
}

// Set parts of Matrix
template <typename T> inline Matrix4<T> &Matrix4<T>::setTranslationPart(const Vector3<T> & translation)	{
	data[12]=translation[0];
	data[13]=translation[1];
	data[14]=translation[2];
	return (*this);
}

template <typename T> Matrix4<T> &Matrix4<T>::setRotationPartEuler(const T angleX, const T angleY, const T angleZ)	{
	T cr = cos( static_cast<T>(basis::PI)*angleX/180.0f );
	T sr = sin( static_cast<T>(basis::PI)*angleX/180.0f );
	T cp = cos( static_cast<T>(basis::PI)*angleY/180.0f );
	T sp = sin( static_cast<T>(basis::PI)*angleY/180.0f );
	T cy = cos( static_cast<T>(basis::PI)*angleZ/180.0f );
	T sy = sin( static_cast<T>(basis::PI)*angleZ/180.0f );

	data[0] = ( cp*cy );
	data[1] = ( cp*sy );
	data[2] = ( -sp );

	T srsp = sr*sp;
	T crsp = cr*sp;

	data[4] = ( srsp*cy-cr*sy );
	data[5] = ( srsp*sy+cr*cy );
	data[6] = ( sr*cp );

	data[8] = ( crsp*cy+sr*sy );
	data[9] = ( crsp*sy-sr*cy );
	data[10] = ( cr*cp );
	return (*this);
}

// Set View Matrix
template <typename T> const Matrix4<T> &Matrix4<T>::setPerspective(T left, T right, T bottom, T top, T n, T f)	{
	
	T nudge=0.999f;		//prevent artifacts with infinite far plane

	zero();

	//check for division by 0
	if(left==right || top==bottom || n==f)
		return (*this);

	data[0]=(2*n)/(right-left);
	data[5]=(2*n)/(top-bottom);

	data[8]=(right+left)/(right-left);
	data[9]=(top+bottom)/(top-bottom);

	if(f != -1)	data[10]=-(f+n)/(f-n);
	else		data[10]=-nudge;		//if f==-1, use an infinite far plane

	data[11]=-1;

	if(f != -1)	data[14]=-(2*f*n)/(f-n);
	else		data[14]=-2*n*nudge;	//if f==-1, use an infinite far plane

	return (*this);
}

template <typename T> const Matrix4<T> &Matrix4<T>::setPerspective(T fovy, T aspect, T n, T f)	{
	T left, right, top, bottom;

	//convert fov from degrees to radians
	fovy *= static_cast<T>(basis::PI)/180.0f;

	top = n * tan(fovy/2.0f);
	bottom = -top;

	left = aspect*bottom;
	right = aspect*top;

	return setPerspective(left, right, bottom, top, n, f);
}

template <typename T> const Matrix4<T> &Matrix4<T>::setOrtho(T left, T right, T bottom, T top, T n, T f)	{
	identity();

	data[0]=2.0f/(right-left);

	data[5]=2.0f/(top-bottom);

	data[10]=-2.0f/(f-n);

	data[12]=-(right+left)/(right-left);
	data[13]=-(top+bottom)/(top-bottom);
	data[14]=-(f+n)/(f-n);

	return (*this);
}

template <typename T> const Matrix4<T> &Matrix4<T>::setView(const Vector3<T> &eye, const Vector3<T> &center, const Vector3<T> &up)	{
	return setView(eye.data, center.data, up.data);
}

template <typename T> const Matrix4<T> &Matrix4<T>::setView(const T eye[3], const T center[3], const T up[3])	{
	identity();

	T x[3], y[3], z[3];

	// Difference eye and center vectors to make Z vector. 
	// Normalize Z. 
	SUB3(z, eye, center);		
	NORMAL3(z);

	// Up vector makes Y vector. 
	CPY3(y, up);
	
	// X vector = Y cross Z. 
	// Recompute Y = Z cross X. 
	// Normalize X,Y 
	CROSS3(x,y,z);	
	CROSS3(y,z,x);	
	NORMAL3(x);
	NORMAL3(y);

	// Build resulting view matrix. 
	data[0+0*4] = x[0];  data[0+1*4] = x[1];	data[0+2*4] = x[2];  
	data[0+3*4] = -x[0]*eye[0] - x[1]*eye[1] - x[2]*eye[2];

	data[1+0*4] = y[0];  data[1+1*4] = y[1];	data[1+2*4] = y[2];  
	data[1+3*4] = -y[0]*eye[0] - y[1]*eye[1] - y[2]*eye[2];

	data[2+0*4] = z[0];  data[2+1*4] = z[1];	data[2+2*4] = z[2];  
	data[2+3*4] = -z[0]*eye[0] - z[1]*eye[1] - z[2]*eye[2];

	data[3+0*4] = 0.0;   data[3+1*4] = 0.0;	data[3+2*4] = 0.0;  data[3+3*4] = 1.0;

	return (*this);
}

//////////////////////////////////////////////////////////////////////////
// Comparison operators
template <typename T> inline const bool Matrix4<T>::operator==(const Matrix4<T>& da) const	{

	for(index_type i=0; i<16; ++i)	{
		if(data[i] != da.data[i])
			return false;
	}
	return true;
}

template <typename T> inline const bool Matrix4<T>::operator!=(const Matrix4<T>& da) const	{
	return !((*this)==da);
}

//////////////////////////////////////////////////////////////////////////
// In place arithmetic
template <typename T> inline Matrix4<T> &Matrix4<T>::operator*=(T d)	{
	return ((*this) = (*this) * d);
}

template <typename T> inline Matrix4<T> &Matrix4<T>::operator/=(T d)	{
	return ((*this) = (*this) / d);
}

template <typename T> inline Matrix4<T> &Matrix4<T>::operator+=(const Matrix4 &da)	{
	return ((*this) = (*this) + da);
}

template <typename T> inline Matrix4<T> &Matrix4<T>::operator-=(const Matrix4 &da)	{
	return ((*this) = (*this) - da);
}

template <typename T> inline Matrix4<T> &Matrix4<T>::operator*=(const Matrix4 &da)	{
	return ((*this) = (*this) * da);
}

//////////////////////////////////////////////////////////////////////////
// Friend methods
template <typename T> inline const Matrix4<T> operator-( const Matrix4<T> &a ) {
	return Matrix4<T>(-a.data[0],-a.data[1],-a.data[2],-a.data[3],
					  -a.data[4],-a.data[5],-a.data[6],-a.data[7],
					  -a.data[8],-a.data[9],-a.data[10],-a.data[11],
					  -a.data[12],-a.data[13],-a.data[14],-a.data[15]);
}
template <typename T> inline const Matrix4<T> operator+( const Matrix4<T> &a, const Matrix4<T> &b ) {
	return Matrix4<T>(a.data[0]+b.data[0],a.data[1]+b.data[1],a.data[2]+b.data[2],a.data[3]+b.data[3],
				  	  a.data[4]+b.data[4],a.data[5]+b.data[5],a.data[6]+b.data[6],a.data[7]+b.data[7],
					  a.data[8]+b.data[8],a.data[9]+b.data[9],a.data[10]+b.data[10],a.data[11]+b.data[11],
					  a.data[12]+b.data[12],a.data[13]+b.data[13],a.data[14]+b.data[14],a.data[15]+b.data[15]);
}
template <typename T> inline const Matrix4<T> operator-( const Matrix4<T> &a, const Matrix4<T> &b ) {
	return Matrix4<T>(a.data[0]-b.data[0],a.data[1]-b.data[1],a.data[2]-b.data[2],a.data[3]-b.data[3],
					  a.data[4]-b.data[4],a.data[5]-b.data[5],a.data[6]-b.data[6],a.data[7]-b.data[7],
					  a.data[8]-b.data[8],a.data[9]-b.data[9],a.data[10]-b.data[10],a.data[11]-b.data[11],
					  a.data[12]-b.data[12],a.data[13]-b.data[13],a.data[14]-b.data[14],a.data[15]-b.data[15]);
}
template <typename T> const Matrix4<T> operator*( const Matrix4<T> &a, const Matrix4<T> &b ) {
	//Optimise for matrices in which bottom row is (0, 0, 0, 1) in both matrices
	if(	a.data[3]==0.0f && a.data[7]==0.0f && a.data[11]==0.0f && a.data[15]==1.0f	&&
		b.data[3]==0.0f && b.data[7]==0.0f && b.data[11]==0.0f && b.data[15]==1.0f)
	{
		return Matrix4<T>(	a.data[0]*b.data[0]+a.data[4]*b.data[1]+a.data[8]*b.data[2],
							a.data[1]*b.data[0]+a.data[5]*b.data[1]+a.data[9]*b.data[2],
							a.data[2]*b.data[0]+a.data[6]*b.data[1]+a.data[10]*b.data[2],
							0.0f,
							a.data[0]*b.data[4]+a.data[4]*b.data[5]+a.data[8]*b.data[6],
							a.data[1]*b.data[4]+a.data[5]*b.data[5]+a.data[9]*b.data[6],
							a.data[2]*b.data[4]+a.data[6]*b.data[5]+a.data[10]*b.data[6],
							0.0f,
							a.data[0]*b.data[8]+a.data[4]*b.data[9]+a.data[8]*b.data[10],
							a.data[1]*b.data[8]+a.data[5]*b.data[9]+a.data[9]*b.data[10],
							a.data[2]*b.data[8]+a.data[6]*b.data[9]+a.data[10]*b.data[10],
							0.0f,
							a.data[0]*b.data[12]+a.data[4]*b.data[13]+a.data[8]*b.data[14]+a.data[12],
							a.data[1]*b.data[12]+a.data[5]*b.data[13]+a.data[9]*b.data[14]+a.data[13],
							a.data[2]*b.data[12]+a.data[6]*b.data[13]+a.data[10]*b.data[14]+a.data[14],
							1.0f);
	}

	//Optimise for when bottom row of 1st matrix is (0, 0, 0, 1)
	if(	a.data[3]==0.0f && a.data[7]==0.0f && a.data[11]==0.0f && a.data[15]==1.0f)
	{
		return Matrix4<T>(	a.data[0]*b.data[0]+a.data[4]*b.data[1]+a.data[8]*b.data[2]+a.data[12]*b.data[3],
							a.data[1]*b.data[0]+a.data[5]*b.data[1]+a.data[9]*b.data[2]+a.data[13]*b.data[3],
							a.data[2]*b.data[0]+a.data[6]*b.data[1]+a.data[10]*b.data[2]+a.data[14]*b.data[3],
							b.data[3],
							a.data[0]*b.data[4]+a.data[4]*b.data[5]+a.data[8]*b.data[6]+a.data[12]*b.data[7],
							a.data[1]*b.data[4]+a.data[5]*b.data[5]+a.data[9]*b.data[6]+a.data[13]*b.data[7],
							a.data[2]*b.data[4]+a.data[6]*b.data[5]+a.data[10]*b.data[6]+a.data[14]*b.data[7],
							b.data[7],
							a.data[0]*b.data[8]+a.data[4]*b.data[9]+a.data[8]*b.data[10]+a.data[12]*b.data[11],
							a.data[1]*b.data[8]+a.data[5]*b.data[9]+a.data[9]*b.data[10]+a.data[13]*b.data[11],
							a.data[2]*b.data[8]+a.data[6]*b.data[9]+a.data[10]*b.data[10]+a.data[14]*b.data[11],
							b.data[11],
							a.data[0]*b.data[12]+a.data[4]*b.data[13]+a.data[8]*b.data[14]+a.data[12]*b.data[15],
							a.data[1]*b.data[12]+a.data[5]*b.data[13]+a.data[9]*b.data[14]+a.data[13]*b.data[15],
							a.data[2]*b.data[12]+a.data[6]*b.data[13]+a.data[10]*b.data[14]+a.data[14]*b.data[15],
							b.data[15]);
	}

	// Optimize for when bottom row of 2nd matrix is (0, 0, 0, 1)
	if(	b.data[3]==0.0f && b.data[7]==0.0f && b.data[11]==0.0f && b.data[15]==1.0f)
	{
		return Matrix4<T>(	a.data[0]*b.data[0]+a.data[4]*b.data[1]+a.data[8]*b.data[2],
							a.data[1]*b.data[0]+a.data[5]*b.data[1]+a.data[9]*b.data[2],
							a.data[2]*b.data[0]+a.data[6]*b.data[1]+a.data[10]*b.data[2],
							a.data[3]*b.data[0]+a.data[7]*b.data[1]+a.data[11]*b.data[2],
							a.data[0]*b.data[4]+a.data[4]*b.data[5]+a.data[8]*b.data[6],
							a.data[1]*b.data[4]+a.data[5]*b.data[5]+a.data[9]*b.data[6],
							a.data[2]*b.data[4]+a.data[6]*b.data[5]+a.data[10]*b.data[6],
							a.data[3]*b.data[4]+a.data[7]*b.data[5]+a.data[11]*b.data[6],
							a.data[0]*b.data[8]+a.data[4]*b.data[9]+a.data[8]*b.data[10],
							a.data[1]*b.data[8]+a.data[5]*b.data[9]+a.data[9]*b.data[10],
							a.data[2]*b.data[8]+a.data[6]*b.data[9]+a.data[10]*b.data[10],
							a.data[3]*b.data[8]+a.data[7]*b.data[9]+a.data[11]*b.data[10],
							a.data[0]*b.data[12]+a.data[4]*b.data[13]+a.data[8]*b.data[14]+a.data[12],
							a.data[1]*b.data[12]+a.data[5]*b.data[13]+a.data[9]*b.data[14]+a.data[13],
							a.data[2]*b.data[12]+a.data[6]*b.data[13]+a.data[10]*b.data[14]+a.data[14],
							a.data[3]*b.data[12]+a.data[7]*b.data[13]+a.data[11]*b.data[14]+a.data[15]);
	}	

	return Matrix4<T>(	a.data[0]*b.data[0]+a.data[4]*b.data[1]+a.data[8]*b.data[2]+a.data[12]*b.data[3],
						a.data[1]*b.data[0]+a.data[5]*b.data[1]+a.data[9]*b.data[2]+a.data[13]*b.data[3],
						a.data[2]*b.data[0]+a.data[6]*b.data[1]+a.data[10]*b.data[2]+a.data[14]*b.data[3],
						a.data[3]*b.data[0]+a.data[7]*b.data[1]+a.data[11]*b.data[2]+a.data[15]*b.data[3],
						a.data[0]*b.data[4]+a.data[4]*b.data[5]+a.data[8]*b.data[6]+a.data[12]*b.data[7],
						a.data[1]*b.data[4]+a.data[5]*b.data[5]+a.data[9]*b.data[6]+a.data[13]*b.data[7],
						a.data[2]*b.data[4]+a.data[6]*b.data[5]+a.data[10]*b.data[6]+a.data[14]*b.data[7],
						a.data[3]*b.data[4]+a.data[7]*b.data[5]+a.data[11]*b.data[6]+a.data[15]*b.data[7],
						a.data[0]*b.data[8]+a.data[4]*b.data[9]+a.data[8]*b.data[10]+a.data[12]*b.data[11],
						a.data[1]*b.data[8]+a.data[5]*b.data[9]+a.data[9]*b.data[10]+a.data[13]*b.data[11],
						a.data[2]*b.data[8]+a.data[6]*b.data[9]+a.data[10]*b.data[10]+a.data[14]*b.data[11],
						a.data[3]*b.data[8]+a.data[7]*b.data[9]+a.data[11]*b.data[10]+a.data[15]*b.data[11],
						a.data[0]*b.data[12]+a.data[4]*b.data[13]+a.data[8]*b.data[14]+a.data[12]*b.data[15],
						a.data[1]*b.data[12]+a.data[5]*b.data[13]+a.data[9]*b.data[14]+a.data[13]*b.data[15],
						a.data[2]*b.data[12]+a.data[6]*b.data[13]+a.data[10]*b.data[14]+a.data[14]*b.data[15],
						a.data[3]*b.data[12]+a.data[7]*b.data[13]+a.data[11]*b.data[14]+a.data[15]*b.data[15]);
}

template <typename T> const Vector4<T> operator*( const Matrix4<T> &a, const Vector4<T> &b ) {
	// Optimize for matrices in which bottom row is (0, 0, 0, 1)

	if(a.data[3]==0.0f && a.data[7]==0.0f && a.data[11]==0.0f && a.data[15]==1.0f)	{
		return Vector4<T>(a.data[0]*b[0] + a.data[4]*b[1] +	a.data[8]*b[2] + a.data[12]*b[3],
						a.data[1]*b[0] + a.data[5]*b[1]	+ a.data[9]*b[2] + a.data[13]*b[3],
						a.data[2]*b[0] + a.data[6]*b[1]	+ a.data[10]*b[2] +	a.data[14]*b[3],
						b[3]);
	}

	return Vector4<T>(	a.data[0]*b[0] + a.data[4]*b[1] + a.data[8]*b[2] + a.data[12]*b[3],
						a.data[1]*b[0] + a.data[5]*b[1]	+ a.data[9]*b[2] + a.data[13]*b[3],
						a.data[2]*b[0] + a.data[6]*b[1] + a.data[10]*b[2]+ a.data[14]*b[3],
						a.data[3]*b[0] + a.data[7]*b[1] + a.data[11]*b[2]+ a.data[15]*b[3]);
}

template <typename T> inline const Vector3<T> operator*( const Matrix4<T> &a, const Vector3<T> &b ) {
	// Optimize for matrices in which bottom row is (0, 0, 0, 1)

	return Vector3<T>(	a.data[0]*b[0] + a.data[4]*b[1] + a.data[8]*b[2] + a.data[12],
						a.data[1]*b[0] + a.data[5]*b[1]	+ a.data[9]*b[2] + a.data[13],
						a.data[2]*b[0] + a.data[6]*b[1] + a.data[10]*b[2]+ a.data[14]);
}

template <typename T> inline const Matrix4<T> operator+( const Matrix4<T> &a, const T b ) {
	return Matrix4<T>( a.data[0]+b,	a.data[1]+b, a.data[2]+b, a.data[3]+b,
					   a.data[4]+b, a.data[5]+b, a.data[6]+b, a.data[7]+b,
					   a.data[8]+b, a.data[9]+b, a.data[10]+b,a.data[11]+b,
					   a.data[12]+b,a.data[13]+b,a.data[14]+b,a.data[15]+b);
}
template <typename T> inline const Matrix4<T> operator-( const Matrix4<T> &a, const T b ) {
	return Matrix4<T>( a.data[0]-b,	a.data[1]-b, a.data[2]-b, a.data[3]-b,
					   a.data[4]-b, a.data[5]-b, a.data[6]-b, a.data[7]-b,
					   a.data[8]-b, a.data[9]-b, a.data[10]-b,a.data[11]-b,
					   a.data[12]-b,a.data[13]-b,a.data[14]-b,a.data[15]-b);
}
template <typename T> inline const Matrix4<T> operator*( const Matrix4<T> &a, const T b ) {
	return Matrix4<T>( a.data[0]*b,	a.data[1]*b, a.data[2]*b, a.data[3]*b,
					   a.data[4]*b, a.data[5]*b, a.data[6]*b, a.data[7]*b,
					   a.data[8]*b, a.data[9]*b, a.data[10]*b,a.data[11]*b,
					   a.data[12]*b,a.data[13]*b,a.data[14]*b,a.data[15]*b);
}
template <typename T> inline const Matrix4<T> operator/( const Matrix4<T> &a, const T b ) {
	return Matrix4<T>( a.data[0]/b,	a.data[1]/b, a.data[2]/b, a.data[3]/b,
					   a.data[4]/b, a.data[5]/b, a.data[6]/b, a.data[7]/b,
					   a.data[8]/b, a.data[9]/b, a.data[10]/b,a.data[11]/b,
					   a.data[12]/b,a.data[13]/b,a.data[14]/b,a.data[15]/b);
}

template <typename T> inline const Matrix4<T> operator+( const T a, const Matrix4<T> &b ) {
	return (b+a);
}
template <typename T> inline const Matrix4<T> operator-( const T a, const Matrix4<T> &b ) {
	return (-b+a);
}
template <typename T> inline const Matrix4<T> operator*( const T a, const Matrix4<T> &b ) {
	return (b*a);
}

template <typename T> std::ostream &operator<<( std::ostream &os, const Matrix4<T> &m ) {
	os	<< "[" << m.data[0] << ", " << m.data[4] << ", " << m.data[8] << ", " << m.data[12] << std::endl
		<< " " << m.data[1] << ", " << m.data[5] << ", " << m.data[9] << ", " << m.data[13] << std::endl
		<< " " << m.data[2] << ", " << m.data[6] << ", " << m.data[10] << ", " << m.data[14] << std::endl
		<< " " << m.data[3] << ", " << m.data[7] << ", " << m.data[11] << ", " << m.data[15] << "]";
	return os;
}

typedef Matrix4<float>			Matrix4f;
typedef Matrix4<double>			Matrix4d;

#undef __ASSERT

}; // namespace arithmetic

#endif	