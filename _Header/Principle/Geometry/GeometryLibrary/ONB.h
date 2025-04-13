#ifndef _GEOMETRY_ONB_H_
#define _GEOMETRY_ONB_H_

#include "../../Arithmetic/Vector/Vector3.h"

namespace geometryspace {

using arithmetic::Vector3;

template<typename T>
class ONB {
public : 
	static const T Epsilon() { return 0.01f; }

public :
	ONB() {}
	ONB(const Vector3<T>& a, const Vector3<T>& b, const Vector3<T>& c) : U(a), V(b), W(c) {}

	//////////////////////////////////////////////////////////////////////////
	// Setter
	void initFromU( const Vector3<T>& u );
	void initFromV( const Vector3<T>& v );
	void initFromW( const Vector3<T>& w );

	void initFromUnitVector( const Vector3<T>& w );

	//////////////////////////////////////////////////////////////////////////
	// Calculate the ONB from two vectors
	// The first one is the Fixed vector (it is just normalized)
	// The second is normalized and its direction can be adjusted
	void initFromUV( const Vector3<T>& u, const Vector3<T>& v );
	void initFromVU( const Vector3<T>& v, const Vector3<T>& u );

	void initFromUW( const Vector3<T>& u, const Vector3<T>& w );
	void initFromWU( const Vector3<T>& w, const Vector3<T>& u );

	void initFromVW( const Vector3<T>& v, const Vector3<T>& w );
	void initFromWV( const Vector3<T>& w, const Vector3<T>& v );

	void setU(const Vector3<T> &v)	{ U = v; }
	void setV(const Vector3<T> &v)	{ V = v; }
	void setW(const Vector3<T> &v)	{ W = v; }
	void set(const Vector3<T>& a, const Vector3<T>& b, const Vector3<T>&c) { U = a; V = b; W = c; }

	//////////////////////////////////////////////////////////////////////////
	// Getter
	const Vector3<T> u() const { return U; }
	const Vector3<T> v() const { return V; }
	const Vector3<T> w() const { return W; }	

	//////////////////////////////////////////////////////////////////////////
	// Friend Method
	template<typename T> friend std::istream &operator>>(std::istream &is, ONB<T>& t);
	template<typename T> friend std::ostream &operator<<(std::ostream &os, const ONB<T>& t);
	template<typename T> friend inline const bool operator==(const ONB<T>& o1, const ONB<T>& o2);	

public :
	Vector3<T> U, V, W;
};

template<typename T> 
void ONB<T>::initFromU( const Vector3<T>& u )	{
	Vector3<T> n(1.0f, 0.0f, 0.0f);
	Vector3<T> m(0.0f, 1.0f, 0.0f);
	U = unitVector(u);
	V = cross(U, n);
	if(arithmetic::sqrmag(V) < ONB<T>::Epsilon())
		V = cross(U,m);

	V.normalize();
	W = cross(U,V);
}

template<typename T> 
void ONB<T>::initFromV( const Vector3<T>& v )	{
	Vector3<T> n(1.0f, 0.0f, 0.0f);
	Vector3<T> m(0.0f, 1.0f, 0.0f);
	V = unitVector(v);
	U = cross(V, n);
	if(arithmetic::sqrmag(U) < ONB<T>::Epsilon())
		U = cross(V,m);

	U.normalize();
	W = cross(U,V);
}

template<typename T> 
void ONB<T>::initFromW( const Vector3<T>& w )	{
	Vector3<T> n(1.0f, 0.0f, 0.0f);
	Vector3<T> m(0.0f, 1.0f, 0.0f);
	W = unitVector(w);
	U = cross(W, n);
	if(arithmetic::sqrmag(U) < ONB<T>::Epsilon())
		U = cross(W,m);

	U.normalize();
	V = cross(W,U);
}

template<typename T> 
void ONB<T>::initFromUnitVector( const Vector3<T>& w )	{
	Vector3<T> n(1.0f, 0.0f, 0.0f);
	Vector3<T> m(0.0f, 1.0f, 0.0f);
	W = w;
	U = cross(W, n);
	if(arithmetic::sqrmag(U) < ONB<T>::Epsilon())
		U = cross(W,m);

	U.normalize();
	V = cross(W,U);
}

template<typename T> 
void ONB<T>::initFromUV( const Vector3<T>& u, const Vector3<T>& v )	{
	U = unitVector(u);
	W = unitVector( cross(u,v) );
	V = cross(W, U);
}

template<typename T> 
void ONB<T>::initFromVU( const Vector3<T>& v, const Vector3<T>& u )	{
	V = unitVector(v);
	W = unitVector( cross(u,v) );
	U = cross(V,W);
}

template<typename T> 
void ONB<T>::initFromUW( const Vector3<T>& u, const Vector3<T>& w )	{
	U = unitVector(u);
	V = unitVector( cross(w,u) );
	W = cross(U,V);
}

template<typename T> 
void ONB<T>::initFromWU( const Vector3<T>& w, const Vector3<T>& u )	{
	W = unitVector(w);
	V = unitVector( cross(w,u) );
	U = cross(V,W);
}

template<typename T> 
void ONB<T>::initFromVW( const Vector3<T>& v, const Vector3<T>& w )	{
	V = unitVector(v);
	U = unitVector( cross(v,w) );
	W = cross(U,V);
}

template<typename T> 
void ONB<T>::initFromWV( const Vector3<T>& w, const Vector3<T>& v )	{
	W = unitVector(w);
	U = unitVector( cross(v,w) );
	V = cross(W,U);
}

template<typename T> 
std::istream& operator>>( std::istream& is, ONB<T>& t )	{
	Vector3<T> new_u, new_v, new_w;
	is >> new_u >> new_v >> new_w;
	t.initFromUV(new_u, new_v);

	return is;
}

template<typename T> 
std::ostream& operator<<( std::ostream& os, const ONB<T>& t ) {
	os << t.U << std::endl << t.V << std::endl << t.W << std::endl;
	return os;
}

template<typename T> 
inline const bool operator==( const ONB<T>& o1, const ONB<T>& o2 ) {
	return (o1.U == o2.U && o1.V == o2.V && o1.W == o2.W);
}

};		// illuminationspace

#endif
