#include <stdafx.h>
#include <Illumination/IlluminationHeader/Core/ONB.h>

using namespace std;

namespace illuminationspace {

void ONB::initFromU( const Vector3<real>& u )	{
	Vector3<real> n(1.0f, 0.0f, 0.0f);
	Vector3<real> m(0.0f, 1.0f, 0.0f);
	U = unitVector(u);
	V = cross(U, n);
	if(arithmetic::sqrmag(V) < ONB::Epsilon())
		V = cross(U,m);

	V.normalize();
	W = cross(U,V);
}

void ONB::initFromV( const Vector3<real>& v )	{
	Vector3<real> n(1.0f, 0.0f, 0.0f);
	Vector3<real> m(0.0f, 1.0f, 0.0f);
	V = unitVector(v);
	U = cross(V, n);
	if(arithmetic::sqrmag(U) < ONB::Epsilon())
		U = cross(V,m);

	U.normalize();
	W = cross(U,V);
}

void ONB::initFromW( const Vector3<real>& w )	{
	Vector3<real> n(1.0f, 0.0f, 0.0f);
	Vector3<real> m(0.0f, 1.0f, 0.0f);
	W = unitVector(w);
	U = cross(W, n);
	if(arithmetic::sqrmag(U) < ONB::Epsilon())
		U = cross(W,m);

	U.normalize();
	V = cross(W,U);
}

void ONB::initFromUnitVector( const Vector3<real>& w )	{
	Vector3<real> n(1.0f, 0.0f, 0.0f);
	Vector3<real> m(0.0f, 1.0f, 0.0f);
	W = w;
	U = cross(W, n);
	if(arithmetic::sqrmag(U) < ONB::Epsilon())
		U = cross(W,m);

	U.normalize();
	V = cross(W,U);
}

void ONB::initFromUV( const Vector3<real>& u, const Vector3<real>& v )	{
	U = unitVector(u);
	W = unitVector( cross(u,v) );
	V = cross(W, U);
}

void ONB::initFromVU( const Vector3<real>& v, const Vector3<real>& u )	{
	V = unitVector(v);
	W = unitVector( cross(u,v) );
	U = cross(V,W);
}

void ONB::initFromUW( const Vector3<real>& u, const Vector3<real>& w )	{
	U = unitVector(u);
	V = unitVector( cross(w,u) );
	W = cross(U,V);
}

void ONB::initFromWU( const Vector3<real>& w, const Vector3<real>& u )	{
	W = unitVector(w);
	V = unitVector( cross(w,u) );
	U = cross(V,W);
}

void ONB::initFromVW( const Vector3<real>& v, const Vector3<real>& w )	{
	V = unitVector(v);
	U = unitVector( cross(v,w) );
	W = cross(U,V);
}

void ONB::initFromWV( const Vector3<real>& w, const Vector3<real>& v )	{
	W = unitVector(w);
	U = unitVector( cross(v,w) );
	V = cross(W,U);
}

std::istream& operator>>( std::istream& is, ONB& t )	{
	Vector3<real> new_u, new_v, new_w;
	is >> new_u >> new_v >> new_w;
	t.initFromUV(new_u, new_v);

	return is;
}

std::ostream& operator<<( std::ostream& os, const ONB& t ) {
	os << t.U << std::endl << t.V << std::endl << t.W << std::endl;
	return os;
}

};	// illuminationspace