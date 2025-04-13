#ifndef _ARITHMETIC_ARRAYUTIL2_H_
#define _ARITHMETIC_ARRAYUTIL2_H_

#include <iostream>
#include "../../../Basis/Basis/basicLibrary.h"

namespace arithmetic {

template<typename T> inline void		PRINT2(std::ostream &os,const T x[2]);
template<typename T> inline void		ZERO2( T x[2]);
template<typename T> inline void		ADD2( T a[2], const T b[2], const T c[2] );
template<typename T> inline void		SUB2( T a[2], const T b[2], const T c[2] );
template<typename T> inline void		AVG2( T a[2], const T b[2], const T c[2] );
template<typename T> inline void		MUL2( T a[2], const T b,    const T c[2] );
template<typename T> inline void		MUL2( T a[2], const T b[2], const T c    );
template<typename T> inline void		NEG2( T a[2], const T b[2] );
template<typename T> inline void		CPY2( T a[2], const T b[2] );
template<typename T> inline const T		SQLENGTH2( const T a[2] );
template<typename T> inline const T		LENGTH2( const T a[2] );
template<typename T> inline const T		SQDIST2( const T a[2], const T b[2] );
template<typename T> inline const T		DIST2( const T a[2], const T b[2] );
template<typename T> void				NORMAL2( T a[2], T eps );
template<typename T> inline const T		DOT2( const T a[2], const T b[2] );
template<typename T> inline void		ROT2( T a[2], const T b[2], const T theta );

template<typename T> inline void		PRINT2x2(std::ostream &os,const T A[2][2]);
template<typename T> inline void		ZERO2x2( T A[2][2] );
template<typename T> inline void		IDENTITY2x2( T A[2][2] );
template<typename T> inline void		ADD2x2( T A[2][2], const T B[2][2], const T C[2][2] );
template<typename T> inline void		SUB2x2( T A[2][2], const T B[2][2], const T C[2][2] );
template<typename T> inline void		NEG2x2( T A[2][2], const T B[2][2] );
template<typename T> inline void		CPY2x2( T A[2][2], const T B[2][2] );
template<typename T> inline void		TRANS2x2( T A[2][2] );
template<typename T> inline void		TRANS2x2( T A[2][2], const T B[2][2] );
template<typename T> inline void		MUL2x2( T R[2][2], const T s      , const T M[2][2] ) ;
template<typename T> inline void		MUL2x2( T R[2][2], const T M[2][2], const T s ) ;
template<typename T> inline void		MUL2x2( T b[2]   , const T A[2][2], const T x[2] );
template<typename T> inline void		MUL2x2( T A[2][2], const T B[2][2], const T C[2][2] ) ;
template<typename T> inline const T		DET2x2( const T a[2], const T b[2] );
template<typename T> inline const T		DET2x2( const T A[2][2] );
template<typename T> const T			INV2x2( T Inv[2][2], const T A[2][2], T eps );
template<typename T> const T			LinSol2x2( T x[2], const T A[2][2], const T b[2], T eps );
template<typename T> inline const T		TRACE2x2( const T A[2][2] );

// 2 VECTOR
template<typename T> inline void PRINT2(std::ostream &os, const T x[2]) {
	os << "[" << x[0] << ", " << x[1] << "]";
	return os;
}

template<typename T> inline void ZERO2( T x[2] )	{
	x[0] = x[1] = 0;
}

template<typename T> inline void ADD2( T a[2], const T b[2], const T c[2] ) { // a = b + c
	a[0] = b[0] + c[0];
	a[1] = b[1] + c[1];
}

inline void ADD2( unsigned int a[2], const unsigned int b[2], const int c[2] ) { // a = b + c
	a[0] = b[0] + c[0];
	a[1] = b[1] + c[1];
}

template<typename T> inline void SUB2( T a[2], const T b[2], const T c[2] ) { // a = b - c
	a[0] = b[0] - c[0];
	a[1] = b[1] - c[1];
}

template<typename T> inline void AVG2( T a[2], const T b[2], const T c[2] ) {
	a[0] = 0.5f * (b[0] + c[0]);
	a[1] = 0.5f * (b[1] + c[1]);
}

template<typename T> inline void MUL2( T a[2], const T b, const T c[2] ) {
	a[0] = b * c[0]; 
	a[1] = b * c[1]; 
}

template<typename T> inline void MUL2( T a[2], const T b[2], const T c    )	{
	a[0] = b[0] * c; 
	a[1] = b[1] * c; 
}

template<typename T> inline void NEG2( T a[2], const T b[2] )	{
	a[0] = -b[0];
	a[1] = -b[1];
}

template<typename T> inline void CPY2( T a[2], const T b[2] )	{
	a[0] = b[0];
	a[1] = b[1];
}

template<typename T> inline const T SQLENGTH2( const T a[2] )	{
	return a[0]*a[0] + a[1]*a[1];
}

template<typename T> inline const T LENGTH2( const T a[2] ) {
	return sqrt(a[0]*a[0] + a[1]*a[1]);
}

template<typename T> inline const T SQDIST2( const T a[2], const T b[2] ) {
	return (a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]);
}

template<typename T> inline const T DIST2( const T a[2], const T b[2] ) {
	return sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]));
}

template<typename T> void NORMAL2( T a[2], T eps = static_cast<T>(EPSILON) ) {

	T length = sqrt( a[0]*a[0]+a[1]*a[1] );

	if( IsAlmostZero( length, eps ) )	{ a[0]=a[1]=0;					}
	else								{ a[0]/=length; a[1]/=length;	}

}

template<typename T> inline const T DOT2( const T a[2], const T b[2] ) {
	return ( a[0]*b[0] + a[1]*b[1] );
}

template<typename T> inline void ROT2( T a[2], const T b[2], const T theta ) {
	a[0] = cos(theta) * b[0]  -  sin(theta) * b[1];
	a[1] = sin(theta) * b[0]  +  cos(theta) * b[1];
}

// 2x2 MATRIX
template<typename T> inline void PRINT2x2(std::ostream &os,const T A[2][2]) {
	os << A[0][0] << " " << A[0][1] << std::endl;
	os << A[1][0] << " " << A[1][1] << std::endl;	
}

template<typename T> inline void ZERO2x2( T A[2][2] ) {
	A[0][0] = A[0][1] = 0;
	A[1][0] = A[1][1] =	0;
}

template<typename T> inline void IDENTITY2x2( T A[2][2] ) {
	A[0][0]=1;	A[0][1]=0;
	A[1][0]=0;	A[1][1]=1;
}

template<typename T> inline void ADD2x2( T A[2][2], const T B[2][2], const T C[2][2] ) {
	A[0][0] = B[0][0] + C[0][0];
	A[0][1] = B[0][1] + C[0][1];
	A[1][0] = B[1][0] + C[1][0];
	A[1][1] = B[1][1] + C[1][1];
}

template<typename T> inline void SUB2x2( T A[2][2], const T B[2][2], const T C[2][2] )  {
	A[0][0] = B[0][0] - C[0][0];
	A[0][1] = B[0][1] - C[0][1];
	A[1][0] = B[1][0] - C[1][0];
	A[1][1] = B[1][1] - C[1][1];
}

template<typename T> inline void NEG2x2( T A[2][2], const T B[2][2] ) {
	A[0][0] = -B[0][0];
	A[0][1] = -B[0][1];
	A[1][0] = -B[1][0];
	A[1][1] = -B[1][1];
}

template<typename T> inline void CPY2x2( T A[2][2], const T B[2][2] ) { // A =  B
	A[0][0] = B[0][0];
	A[0][1] = B[0][1];
	A[1][0] = B[1][0];
	A[1][1] = B[1][1];
}

template<typename T> inline void TRANS2x2( T A[2][2] ) { // A = A^T
	T temp = A[0][1];	A[0][1] = A[1][0];	A[1][0] = temp;
}

template<typename T> inline void TRANS2x2( T A[2][2], const T B[2][2] ) { // A = B^T
	A[0][0] = B[0][0];
	A[0][1] = B[1][0];
	A[1][0] = B[0][1];
	A[1][1] = B[1][1];
}

template<typename T> inline void MUL2x2( T R[2][2], const T s, const T M[2][2] ) {
	R[0][0] = s*M[0][0];
	R[0][1] = s*M[0][1];
	R[1][0] = s*M[1][0];
	R[1][1] = s*M[1][1];
}

template<typename T> inline void MUL2x2( T R[2][2], const T M[2][2] , const T s ) {
	MUL2x2(R, s, M);
}

template<typename T> inline void MUL2x2( T b[2], const T A[2][2], const T x[2] ) { // b = Ax
	b[0] = A[0][0]*x[0] + A[0][1]*x[1];
	b[1] = A[1][0]*x[0] + A[1][1]*x[1];
}

template<typename T> inline void MUL2x2( T A[2][2], const T B[2][2], const T C[2][2] ) { // A = BC
	A[0][0] = B[0][0]*C[0][0] + B[0][1]*C[1][0];
	A[0][1] = B[0][0]*C[0][1] + B[0][1]*C[1][1];
	A[1][0] = B[1][0]*C[0][0] + B[1][1]*C[1][0];
	A[1][1] = B[1][0]*C[0][1] + B[1][1]*C[1][1];
}

template<typename T> inline void MUL2x2( T** A, const T B[2][2], const T C[2][2] ) { // A = BC
	A[0][0] = B[0][0]*C[0][0] + B[0][1]*C[1][0];
	A[0][1] = B[0][0]*C[0][1] + B[0][1]*C[1][1];
	A[1][0] = B[1][0]*C[0][0] + B[1][1]*C[1][0];
	A[1][1] = B[1][0]*C[0][1] + B[1][1]*C[1][1];
}

template<typename T> inline const T DET2x2( T a[2], const T b[2] ) { 	
	return a[0]*b[1] - a[1]*b[0];
}

template<typename T> inline const T DET2x2( const T A[2][2] ) {
	return A[0][0]*A[1][1] - A[0][1]*A[1][0];
}

template<typename T> const T INV2x2( T Inv[2][2], const T A[2][2], T eps = static_cast<T>(EPSILON) ) { // Inv = A^(-1)

	T det = DET2x2(A);

	if( IsAlmostZero(det, eps) ) { // if not exists
		ZERO2x2(Inv);
		return 0;
	}

	T invDet = 1 / det;

	Inv[0][0] =  A[1][1] * invDet;
	Inv[0][1] = -A[0][1] * invDet;
	Inv[1][0] = -A[1][0] * invDet;
	Inv[1][1] =  A[0][0] * invDet;

	return det;

}

template<typename T> const T LinSol2x2( T x[2], const T A[2][2], const T b[2], T eps = static_cast<T>(EPSILON) ) {

	T det = A[0][0]*A[1][1]-A[0][1]*A[1][0];

	if(IsAlmostZero(det, eps)) { x[0] = x[1] = 0; return det; }

	x[0] = ( A[1][1] * b[0]  -  A[0][1] * b[1] ) / det;
	x[1] = ( A[0][0] * b[1]  -  A[1][0] * b[0] ) / det;

	return det;

}

template<typename T> inline const T TRACE2x2( const T A[2][2] )	{
	return A[0][0] + A[1][1];
}

}; // namespace arithmetic

#endif