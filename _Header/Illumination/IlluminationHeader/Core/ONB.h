#ifndef _ILLUMINATION_ONB_H_
#define _ILLUMINATION_ONB_H_

#include "../../../Principle/Arithmetic/Vector/Vector3.h"
#include "../../IlluminationSpace.h"
#include "../Illumination/Environments.h"

namespace illuminationspace {

using arithmetic::Vector3;

class ILLUMINATIONDLL ONB {
public : 
	static const real Epsilon() { return 0.01f; }

public :
	ONB() {}
	ONB(const Vector3<real>& a, const Vector3<real>& b, const Vector3<real>& c) : U(a), V(b), W(c) {}

	//////////////////////////////////////////////////////////////////////////
	// Setter
	void initFromU( const Vector3<real>& u );
	void initFromV( const Vector3<real>& v );
	void initFromW( const Vector3<real>& w );

	void initFromUnitVector( const Vector3<real>& w );

	//////////////////////////////////////////////////////////////////////////
	// Calculate the ONB from two vectors
	// The first one is the Fixed vector (it is just normalized)
	// The second is normalized and its direction can be adjusted
	void initFromUV( const Vector3<real>& u, const Vector3<real>& v );
	void initFromVU( const Vector3<real>& v, const Vector3<real>& u );

	void initFromUW( const Vector3<real>& u, const Vector3<real>& w );
	void initFromWU( const Vector3<real>& w, const Vector3<real>& u );

	void initFromVW( const Vector3<real>& v, const Vector3<real>& w );
	void initFromWV( const Vector3<real>& w, const Vector3<real>& v );

	void setU(const Vector3<real> &v)	{ U = v; }
	void setV(const Vector3<real> &v)	{ V = v; }
	void setW(const Vector3<real> &v)	{ W = v; }
	void set(const Vector3<real>& a, const Vector3<real>& b, const Vector3<real>&c) { U = a; V = b; W = c; }

	//////////////////////////////////////////////////////////////////////////
	// Getter
	const Vector3<real> u() const { return U; }
	const Vector3<real> v() const { return V; }
	const Vector3<real> w() const { return W; }	

	//////////////////////////////////////////////////////////////////////////
	// Friend Method
	friend ILLUMINATIONDLL std::istream &operator>>(std::istream &is, ONB &t);
	friend ILLUMINATIONDLL std::ostream &operator<<(std::ostream &os, const ONB &t);
	friend inline const bool operator==(const ONB& o1, const ONB &o2);	

public :
	Vector3<real> U, V, W;
};

inline const bool operator==( const ONB& o1, const ONB& o2 ) {
	return (o1.U == o2.U && o1.V == o2.V && o1.W == o2.W);
}

};		// illuminationspace

#endif
