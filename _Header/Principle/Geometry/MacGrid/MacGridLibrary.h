#ifndef _GEOMETRYSPACE_MACGRID_LIBRARY_H_
#define _GEOMETRYSPACE_MACGRID_LIBRARY_H_

#include "../../Arithmetic/GridLibrary/Grid2.h"
#include "../../Arithmetic/GridLibrary/Grid3.h"

namespace geometryspace {

using arithmetic::Grid2;
using arithmetic::Grid3;

template<typename T>
class MacGridLibrary {
public : 
	typedef typename Grid2<T,T>::index_type index_type;

	//////////////////////////////////////////////////////////////////////////
	// Phi Manipulation
	static void getNormals(Grid2<Vector2<T>,T>&, const Grid2<T,T> &phi);
	static void getNormals(Grid3<Vector3<T>,T>&, const Grid3<T,T> &phi);
	static void getCurvatures(Grid2<T,T>&, const Grid2<T,T> &phi);
	static void getCurvatures(Grid3<T,T>&, const Grid3<T,T> &phi);

	static const Vector2<T> getNormal(const Grid2<T,T> &phi, index_type, index_type);
	static const Vector3<T> getNormal(const Grid3<T,T> &phi, index_type, index_type, index_type);
	static const T getCurvature(const Grid2<T,T> &phi, index_type, index_type);
	static const T getCurvature(const Grid3<T,T> &phi, index_type, index_type, index_type);

	static const T delta( T _phi, T epsilon = 1.5f );
	static const T heaviside( T _phi, T _valPos, T _valNeg, T epsilon = 1.5f );
};

template<typename T>
void MacGridLibrary<T>::getNormals(Grid2<Vector2<T>, T> &normals, const Grid2<T,T> &phi) {
	for( index_type i = 1; i < phi.nX-1; ++i ) 
		for ( index_type  j = 1; j < phi.nY-1; ++j ) 
			normals(i,j) = getNormal(phi, i, j);

	normals.setNeumannBoundary();	
}

template<typename T>
void MacGridLibrary<T>::getNormals(Grid3<Vector3<T>, T> &normals, const Grid3<T,T> &phi) {

#ifdef USE_OPENMP
	int my_begin = int(1);
	int my_end = int(phi.nX-1);
	int i;

#pragma omp parallel for
#else
	index_type my_begin = 1;
	index_type my_end = phi.nX-1;
	index_type i;
#endif
	for( i = my_begin; i < my_end; ++i ) 
		for ( index_type  j = 1; j < phi.nY-1; ++j ) 
			for ( index_type  k = 1; k < phi.nZ-1; ++k ) 
				normals(i,j,k) = getNormal(phi, i, j, k);

	normals.setNeumannBoundary();	
}

template<typename T>
void MacGridLibrary<T>::getCurvatures(Grid2<T,T> &curvatures, const Grid2<T,T> &phi) {
	for( index_type i = 1; i < phi.nX-1; ++i ) 
		for ( index_type j = 1; j < phi.nY-1; ++j ) 
			curvatures(i,j) = getCurvature(phi, i, j);

	curvatures.setNeumannBoundary();
}

template<typename T>
void MacGridLibrary<T>::getCurvatures(Grid3<T,T> &curvature, const Grid3<T,T> &phi) {

#ifdef USE_OPENMP
	int my_begin = int(1);
	int my_end = int(phi.nX-1);
	int i;

#pragma omp parallel for
#else
	typename Grid3<T,T>::index_type my_begin = 1;
	typename Grid3<T,T>::index_type my_end = phi.nX-1;
	typename Grid3<T,T>::index_type i;
#endif
	for( i = my_begin; i < my_end; ++i ) 
		for(typename Grid3<T,T>::index_type j=1;j<phi.nY-1;++j) 
			for(typename Grid3<T,T>::index_type k=1;k<phi.nZ-1;++k) 
				curvature(i,j,k) = getCurvature(phi, i, j, k);

	curvature.setNeumannBoundary();
}

template<typename T>
const Vector2<T> MacGridLibrary<T>::getNormal(const Grid2<T,T> &phi, index_type i, index_type j) {
	Vector2<T> ret(phi.dxInvHalf * (phi(i+1,j)-phi(i-1,j))
				 , phi.dyInvHalf * (phi(i,j+1)-phi(i,j-1)) );

	arithmetic::normalize(ret);

	return ret;
}

template<typename T>
const Vector3<T> MacGridLibrary<T>::getNormal(const Grid3<T,T> &phi, index_type i, index_type j, index_type k) {
	Vector3<T> ret(phi.dxInvHalf * (phi(i+1,j,k)-phi(i-1,j,k))
				 , phi.dyInvHalf * (phi(i,j+1,k)-phi(i,j-1,k))
				 , phi.dzInvHalf * (phi(i,j,k+1)-phi(i,j,k-1)) );
	
	arithmetic::normalize(ret);

	return ret;
}

template<typename T>
const T MacGridLibrary<T>::getCurvature(const Grid2<T,T> &phi, index_type i, index_type j) {
	static const T eps = 1.0e-8;
	static const T two = T(2.0);

	const T hInv = phi.dxInv, hInvHalf = phi.dxInvHalf;
	const T hInvSQ = hInv * hInv;
	const T hInvHalfSQ = hInvHalf * hInvHalf;

	T c = 0;
	T phix = hInvHalf * (phi(i+1,j) - phi(i-1,j));
	T phiy = hInvHalf * (phi(i,j+1) - phi(i,j-1));

	T phixx = hInvSQ * (phi(i+1,j) - two * phi(i,j) + phi(i-1,j));
	T phiyy = hInvSQ * (phi(i,j+1) - two * phi(i,j) + phi(i,j-1));

	T phixy = hInvHalfSQ * (phi(i+1,j+1) - phi(i+1,j-1) - phi(i-1,j+1) + phi(i-1,j-1));

	if ( ::abs(phix*phix+phiy*phiy) > eps )
		c = (phix*phix*phiyy - two*phix*phiy*phixy + phiy*phiy*phixx) / CUBE(sqrt(phix*phix+phiy*phiy));

	if ( ::abs(c) > hInv ) 
		c = SIGN(c) * hInv;

	return c;
}

template<typename T>
const T MacGridLibrary<T>::getCurvature(const Grid3<T,T> &phi, index_type i, index_type j, index_type k) {
	static const T eps = 1.0e-8;
	static const T two = T(2.0);

	const T hInv = phi.dxInv, hInvHalf = phi.dxInvHalf;
	const T hInvSQ = hInv * hInv;
	const T hInvHalfSQ = hInvHalf * hInvHalf;

	T c = 0;
	const T phix = hInvHalf * (phi(i+1,j,k)-phi(i-1,j,k));
	const T phiy = hInvHalf * (phi(i,j+1,k)-phi(i,j-1,k));
	const T phiz = hInvHalf * (phi(i,j,k+1)-phi(i,j,k-1));

	const T phixy = hInvHalfSQ * (phi(i+1,j+1,k)-phi(i+1,j-1,k)-phi(i-1,j+1,k)+phi(i-1,j-1,k));
	const T phixz = hInvHalfSQ * (phi(i+1,j,k+1)-phi(i+1,j,k-1)-phi(i-1,j,k+1)+phi(i-1,j,k-1));
	const T phiyz = hInvHalfSQ * (phi(i,j+1,k+1)-phi(i,j+1,k-1)-phi(i,j-1,k+1)+phi(i,j-1,k-1));

	const T phixx = hInvSQ * (phi(i+1,j,k) - two * phi(i,j,k) + phi(i-1,j,k));
	const T phiyy = hInvSQ * (phi(i,j+1,k) - two * phi(i,j,k) + phi(i,j-1,k));
	const T phizz = hInvSQ * (phi(i,j,k+1) - two * phi(i,j,k) + phi(i,j,k-1));

	if ( ::abs(phix*phix+phiy*phiy+phiz*phiz) > eps )
		c = ( SQ(phix)*phiyy - two*phix*phiy*phixy + SQ(phiy)*phixx + SQ(phix)*phizz - two*phix*phiz*phixz + SQ(phiz)*phixx	
			+ SQ(phiy)*phizz - two*phiy*phiz*phiyz + SQ(phiz)*phiyy ) 
			/ CUBE(sqrt(phix*phix+phiy*phiy+phiz*phiz));

	if ( ::abs(c) > hInv ) 
		c = SIGN(c) * hInv;

	return c;
}

template<typename T>
const T MacGridLibrary<T>::delta( T _phi, T epsilon ) {
	if ( _phi > epsilon || _phi < -epsilon )
		return T(0.0);		
	else
		return T(0.5)/epsilon * ( T(1.0) + cos( PI*_phi/epsilon ) );		
}

template<typename T>
const T MacGridLibrary<T>::heaviside( T _phi, T _valPos, T _valNeg, T epsilon ) {
	if ( _phi > epsilon ) {
		return _valPos;
	} else {
		if ( _phi < -epsilon )
			return _valNeg;
		else
			return (_valNeg+(_valPos - _valNeg)*(T(0.5)+T(0.5)*_phi/(epsilon)+T(0.5)/PI*sin(PI*_phi/(epsilon))));
	}
}

};		// geometryspace

#endif