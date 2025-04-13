#ifndef _ARITHMETIC_ARRAYUTIL3_H_
#define _ARITHMETIC_ARRAYUTIL3_H_

#include <iostream>
#include "../../../Basis/Basis/basicLibrary.h"

namespace arithmetic {

template<typename T> inline void		PRINT3(std::ostream &os,const T x[3]);
template<typename T> inline void		ZERO3( T x[3]);
template<typename T> inline void		ADD3( T a[3], const T b[3], const T c[3] );
template<typename T> inline void		SUB3( T a[3], const T b[3], const T c[3] );
template<typename T> inline void		AVG3( T a[3], const T b[3], const T c[3] );
template<typename T> inline void		MUL3( T a[3], const T b,    const T c[3] );
template<typename T> inline void		MUL3( T a[3], const T b[3], const T c    );
template<typename T> inline void		NEG3( T a[3], const T b[3] );
template<typename T> inline void		CPY3( T a[3], const T b[3] );
template<typename T> inline const T		SQLENGTH3( const T a[3] );
template<typename T> inline const T		LENGTH3( const T a[3] );
template<typename T> inline const T		DIST3( const T a[3], const T b[3] );
template<typename T> inline const T		SQDIST3( const T a[3], const T b[3] );
template<typename T> void				NORMAL3( T a[3], T eps );
template<typename T> const T			ANGLE3( const T a[3], const T b[3]);
template<typename T> inline const T		DOT3( const T a[3], const T b[3] );
template<typename T> inline void		CROSS3( T a[3], const T b[3], const T c[3] );

template<typename T> inline void		PRINT3(std::ostream &os,const T A[3][3]);
template<typename T> inline void		ZERO3x3( T A[3][3] );
template<typename T> inline void		IDENTITY3x3( T A[3][3] );
template<typename T> inline void		ADD3x3( T A[3][3], const T B[3][3], const T C[3][3] );
template<typename T> inline void		SUB3x3( T A[3][3], const T B[3][3], const T C[3][3] );
template<typename T> inline void		NEG3x3( T A[3][3], const T B[3][3] );
template<typename T> inline void		CPY3x3( T A[3][3], const T B[3][3] );
template<typename T> inline void		TRANS3x3( T A[3][3] );
template<typename T> inline void		TRANS3x3( T A[3][3], const T B[3][3] );
template<typename T> inline void		MUL3x3( T R[3][3], const T s      , const T M[3][3] );
template<typename T> inline void		MUL3x3( T R[3][3], const T M[3][3], const T s       );
template<typename T> inline void		MUL3x3( T b[3]   , const T A[3][3], const T x[3]	 );
template<typename T> inline void		MUL3x3( T A[3][3], const T B[3][3], const T C[3][3] );
template<typename T> inline const T		DET3x3( const T a[3], const T b[3], const T c[3] );
template<typename T> inline const T		DET3x3( const T A[3][3] );
template<typename T> const T			INV3x3( T Inv[3][3], const T A[3][3], T eps );
template<typename T> const T			LinSol3x3( T x[3], const T A[3][3], const T b[3], T eps );
template<typename T> inline const T		TRACE3x3( const T A[3][3] );

// 3 VECTOR
template<typename T> inline void PRINT3(std::ostream &os,const T x[3]) {
	os << "[" << x[0] << ", " << x[1] << ", " << x[2] << "]";
	return os;
}
template<typename T> inline void ZERO3( T x[3] )	{
	x[0] = x[1] = x[2] = 0;
}

template<typename T> inline void ADD3( T a[3], const T b[3], const T c[3] ) { // a = b + c
	a[0] = b[0] + c[0];
	a[1] = b[1] + c[1];
	a[2] = b[2] + c[2];
}

inline void ADD3( unsigned int a[3], const unsigned int b[3], const int c[3] ) { // a = b + c
	a[0] = b[0] + c[0];
	a[1] = b[1] + c[1];
	a[2] = b[2] + c[2];
}

template<typename T> inline void SUB3( T a[3], const T b[3], const T c[3] ) { // a = b - c
	a[0] = b[0] - c[0];
	a[1] = b[1] - c[1];
	a[2] = b[2] - c[2];
}

template<typename T> inline void AVG3( T a[3], const T b[3], const T c[3] ) {
	a[0] = 0.5f * (b[0] + c[0]);
	a[1] = 0.5f * (b[1] + c[1]);
	a[2] = 0.5f * (b[2] + c[2]);
}

template<typename T> inline void MUL3( T a[3], const T b,    const T c[3] ) {
	a[0] = b * c[0]; 
	a[1] = b * c[1]; 
	a[2] = b * c[2];
}

template<typename T> inline void MUL3( T a[3], const T b[3], const T c    )	{
	a[0] = b[0] * c; 
	a[1] = b[1] * c; 
	a[2] = b[2] * c;
}

template<typename T> inline void NEG3( T a[3], const T b[3] )	{
	a[0] = -b[0];
	a[1] = -b[1];
	a[2] = -b[2];
}

template<typename T> inline void CPY3( T a[3], const T b[3] )	{
	a[0] = b[0];
	a[1] = b[1];
	a[2] = b[2];
}

template<typename T> inline const T SQLENGTH3( const T a[3] ) {
	return a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
}

template<typename T> inline const T LENGTH3( const T a[3] ) {
	return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

template<typename T> inline const T SQDIST3( const T a[3], const T b[3] ) {
	return (a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]) ;
}

template<typename T> inline const T DIST3( const T a[3], const T b[3] ) {
	return sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]));
}

template<typename T> void NORMAL3( T a[3], T eps = static_cast<T>(EPSILON) ) {

	T length = sqrt( a[0]*a[0]+a[1]*a[1]+a[2]*a[2] );

	if( IsAlmostZero( length, eps ) )	{ a[0]=a[1]=a[2]=0;                         }
	else								{ a[0]/=length; a[1]/=length; a[2]/=length; }

}

template<typename T> const T ANGLE3( const T a[3], const T b[3] ) {

	T dotProduct = DOT3(a,b);
	T AB = sqrt(DOT3(a,a)*DOT3(b,b));
	T angle = acos(dotProduct/AB);

	if(_isnan(angle)) {
		std::cout << "Dangerous Angle " << angle << "\n";
		std::cout << dotProduct/AB << "\n";
		return 0;
	}

	return angle;
}

template<typename T> inline const T DOT3( const T a[3], const T b[3] ) {
	return ( a[0]*b[0] + a[1]*b[1] + a[2]*b[2] );
}

template<typename T> inline void CROSS3( T a[3], const T b[3], const T c[3] ) {
	a[0] = b[1]*c[2] - b[2]*c[1];
	a[1] = b[2]*c[0] - b[0]*c[2];
	a[2] = b[0]*c[1] - b[1]*c[0];

}

// 3x3 MATRIX
template<typename T> inline void PRINT3x3(std::ostream &os,const T A[3][3]) {
	os << A[0][0] << " " << A[0][1] << " " << A[0][2] << std::endl;
	os << A[1][0] << " " << A[1][1] << " " << A[1][2] << std::endl;
	os << A[2][0] << " " << A[2][1] << " " << A[2][2] << std::endl;
}

template<typename T> inline void ZERO3x3( T A[3][3] ) {
	A[0][0] = A[0][1] = A[0][2] = 0;
	A[1][0] = A[1][1] =	A[1][2] = 0;
	A[2][0] = A[2][1] =	A[2][2] = 0;
}

template<typename T> inline void IDENTITY3x3( T A[3][3] ) {
	A[0][0]=1;	A[0][1]=0;	A[0][2]=0;
	A[1][0]=0;	A[1][1]=1;	A[1][2]=0;
	A[2][0]=0;	A[2][1]=0;	A[2][2]=1;
}

template<typename T> inline void ADD3x3( T A[3][3], const T B[3][3], const T C[3][3] ) {
	A[0][0] = B[0][0] + C[0][0];
	A[0][1] = B[0][1] + C[0][1];
	A[0][2] = B[0][2] + C[0][2];
	A[1][0] = B[1][0] + C[1][0];
	A[1][1] = B[1][1] + C[1][1];
	A[1][2] = B[1][2] + C[1][2];
	A[2][0] = B[2][0] + C[2][0];
	A[2][1] = B[2][1] + C[2][1];
	A[2][2] = B[2][2] + C[2][2];
}

template<typename T> inline void SUB3x3( T A[3][3], const T B[3][3], const T C[3][3] )  {
	A[0][0] = B[0][0] - C[0][0];
	A[0][1] = B[0][1] - C[0][1];
	A[0][2] = B[0][2] - C[0][2];
	A[1][0] = B[1][0] - C[1][0];
	A[1][1] = B[1][1] - C[1][1];
	A[1][2] = B[1][2] - C[1][2];
	A[2][0] = B[2][0] - C[2][0];
	A[2][1] = B[2][1] - C[2][1];
	A[2][2] = B[2][2] - C[2][2];
}

template<typename T> inline void NEG3x3( T A[3][3], const T B[3][3] ) {
	A[0][0] = -B[0][0];
	A[0][1] = -B[0][1];
	A[0][2] = -B[0][2];
	A[1][0] = -B[1][0];
	A[1][1] = -B[1][1];
	A[1][2] = -B[1][2];
	A[2][0] = -B[2][0];
	A[2][1] = -B[2][1];
	A[2][2] = -B[2][2];
}

template<typename T> inline void CPY3x3( T A[3][3], const T B[3][3] ) { // A =  B

	A[0][0] = B[0][0];
	A[0][1] = B[0][1];
	A[0][2] = B[0][2];

	A[1][0] = B[1][0];
	A[1][1] = B[1][1];
	A[1][2] = B[1][2];

	A[2][0] = B[2][0];
	A[2][1] = B[2][1];
	A[2][2] = B[2][2];

}

template<typename T> inline void TRANS3x3( T A[3][3] ) { // A = A^T

	T temp;

	temp = A[0][1];	A[0][1] = A[1][0];	A[1][0] = temp;
	temp = A[0][2];	A[0][2] = A[2][0];	A[2][0] = temp;
	temp = A[1][2];	A[1][2] = A[2][1];	A[2][1] = temp;

}

template<typename T> inline void TRANS3x3( T A[3][3], const T B[3][3] ) { // A = B^T

	A[0][0] = B[0][0];
	A[0][1] = B[1][0];
	A[0][2] = B[2][0];
	A[1][0] = B[0][1];
	A[1][1] = B[1][1];
	A[1][2] = B[2][1];
	A[2][0] = B[0][2];
	A[2][1] = B[1][2];
	A[2][2] = B[2][2];

}

template<typename T> inline void MUL3x3( T R[3][3], const T s, const T M[3][3] ) {

	R[0][0] = s*M[0][0];
	R[0][1] = s*M[0][1];
	R[0][2] = s*M[0][2];

	R[1][0] = s*M[1][0];
	R[1][1] = s*M[1][1];
	R[1][2] = s*M[1][2];

	R[2][0] = s*M[2][0];
	R[2][1] = s*M[2][1];
	R[2][2] = s*M[2][2];

}

template<typename T> inline void MUL3x3( T R[3][3], const T M[3][3] , const T s ) {
	MUL3x3(R, s, M);
}

template<typename T> inline void MUL3x3( T b[3], const T A[3][3], const T x[3] ) { // b = Ax

	b[0] = A[0][0]*x[0] + A[0][1]*x[1] + A[0][2]*x[2];
	b[1] = A[1][0]*x[0] + A[1][1]*x[1] + A[1][2]*x[2];
	b[2] = A[2][0]*x[0] + A[2][1]*x[1] + A[2][2]*x[2];

}

template<typename T> inline void MUL3x3( T A[3][3], const T B[3][3], const T C[3][3] ) { // A = BC

	A[0][0] = B[0][0]*C[0][0] + B[0][1]*C[1][0] + B[0][2]*C[2][0];
	A[0][1] = B[0][0]*C[0][1] + B[0][1]*C[1][1] + B[0][2]*C[2][1];
	A[0][2] = B[0][0]*C[0][2] + B[0][1]*C[1][2] + B[0][2]*C[2][2];

	A[1][0] = B[1][0]*C[0][0] + B[1][1]*C[1][0] + B[1][2]*C[2][0];
	A[1][1] = B[1][0]*C[0][1] + B[1][1]*C[1][1] + B[1][2]*C[2][1];
	A[1][2] = B[1][0]*C[0][2] + B[1][1]*C[1][2] + B[1][2]*C[2][2];

	A[2][0] = B[2][0]*C[0][0] + B[2][1]*C[1][0] + B[2][2]*C[2][0];
	A[2][1] = B[2][0]*C[0][1] + B[2][1]*C[1][1] + B[2][2]*C[2][1];
	A[2][2] = B[2][0]*C[0][2] + B[2][1]*C[1][2] + B[2][2]*C[2][2];

}

template<typename T> inline void MUL3x3( T** A, const T B[3][3], const T C[3][3] ) { // A = BC

	A[0][0] = B[0][0]*C[0][0] + B[0][1]*C[1][0] + B[0][2]*C[2][0];
	A[0][1] = B[0][0]*C[0][1] + B[0][1]*C[1][1] + B[0][2]*C[2][1];
	A[0][2] = B[0][0]*C[0][2] + B[0][1]*C[1][2] + B[0][2]*C[2][2];

	A[1][0] = B[1][0]*C[0][0] + B[1][1]*C[1][0] + B[1][2]*C[2][0];
	A[1][1] = B[1][0]*C[0][1] + B[1][1]*C[1][1] + B[1][2]*C[2][1];
	A[1][2] = B[1][0]*C[0][2] + B[1][1]*C[1][2] + B[1][2]*C[2][2];

	A[2][0] = B[2][0]*C[0][0] + B[2][1]*C[1][0] + B[2][2]*C[2][0];
	A[2][1] = B[2][0]*C[0][1] + B[2][1]*C[1][1] + B[2][2]*C[2][1];
	A[2][2] = B[2][0]*C[0][2] + B[2][1]*C[1][2] + B[2][2]*C[2][2];

}

template<typename T> inline const T DET3x3( const T a[3], const T b[3], const T c[3] ) { // Get the volume of the parallelepiped spanned by three independent vectors.

	return ( a[0]*( b[1]*c[2] - c[1]*b[2] )
		- b[0]*( a[1]*c[2] - c[1]*a[2] )
		+ c[0]*( a[1]*b[2] - b[1]*a[2] ) );

}

template<typename T> inline const T DET3x3( const T A[3][3] ) {

	return ( A[0][0]*( A[1][1]*A[2][2] - A[1][2]*A[2][1] )
		- A[0][1]*( A[1][0]*A[2][2] - A[1][2]*A[2][0] )
		+ A[0][2]*( A[1][0]*A[2][1] - A[1][1]*A[2][0] ) );

}

template<typename T> const T INV3x3( T Inv[3][3], const T A[3][3], T eps = static_cast<T>(EPSILON) ) { // Inv = A^(-1)

	T det = DET3x3(A);

	if( IsAlmostZero(det, eps) ) { // if not exists
		ZERO3x3(Inv);
		return 0;
	}

	T invDET = 1 / det;

	Inv[0][0] = (  A[1][1] * A[2][2]  -  A[1][2] * A[2][1] ) * invDET;
	Inv[0][1] = ( -A[0][1] * A[2][2]  +  A[0][2] * A[2][1] ) * invDET;
	Inv[0][2] = (  A[0][1] * A[1][2]  -  A[0][2] * A[1][1] ) * invDET;
	Inv[1][0] = ( -A[1][0] * A[2][2]  +  A[1][2] * A[2][0] ) * invDET;
	Inv[1][1] = (  A[0][0] * A[2][2]  -  A[0][2] * A[2][0] ) * invDET;
	Inv[1][2] = ( -A[0][0] * A[1][2]  +  A[0][2] * A[1][0] ) * invDET;
	Inv[2][0] = (  A[1][0] * A[2][1]  -  A[1][1] * A[2][0] ) * invDET;
	Inv[2][1] = ( -A[0][0] * A[2][1]  +  A[0][1] * A[2][0] ) * invDET;
	Inv[2][2] = (  A[0][0] * A[1][1]  -  A[0][1] * A[1][0] ) * invDET;

	return det;

}

template<typename T> const T LinSol3x3( T x[3], const T A[3][3], const T b[3], T eps = static_cast<T>(EPSILON) ) {

	T d1, d2, d3, d;

	T B[3][3];

	TRANS3x3(B, A);
	d  = DET3x3(B[0], B[1], B[2]);

	if( IsAlmostZero(d, eps) ) { x[0]=x[1]=x[2]=0; return d; }

	d1 = DET3x3( b,    B[1], B[2] );
	d2 = DET3x3( B[0], b,    B[2] );
	d3 = DET3x3( B[0], B[1], b    );

	x[0] = d1/d;
	x[1] = d2/d;
	x[2] = d3/d;

	return d;

}

template<typename T> inline const T TRACE3x3( const T A[3][3] )	{
	return A[0][0] + A[1][1] + A[2][2];
}

};	 // namespace arithmetic

#endif