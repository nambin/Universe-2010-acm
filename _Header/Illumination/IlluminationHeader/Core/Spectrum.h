#ifndef _ILLUMINATION_SPECTRUM_H_
#define _ILLUMINATION_SPECTRUM_H_

#include <cmath>
#include <iostream>
#include "../../../Basis/Color/RGBA.h"
#include "../../IlluminationSpace.h"
#include "../Illumination/Constants.h"

namespace illuminationspace {

#ifdef _DEBUG
#define __ASSERT(c) assert(c)
#else
#define __ASSERT(c) {}
#endif

using std::ostream;
using namespace basis;

#define N_COLOR_SAMPLES 3

class ILLUMINATIONDLL Spectrum	{
public :	
	//////////////////////////////////////////////////////////////////////////
	// Constructors
	inline explicit Spectrum(real v = 0.0f);
	inline explicit Spectrum(real cs[N_COLOR_SAMPLES]);
	inline Spectrum(real c0, real c1, real c2);
	inline Spectrum(const RGBA &rgba);

	void init() { set(0.0f); }

	//////////////////////////////////////////////////////////////////////////
	// Assignment and set	
	inline void set(real v);
	inline void set(real cs[N_COLOR_SAMPLES]);
	inline void set(real c0, real c1, real c2);
	inline void set(const Spectrum &rhs);
	inline void set(const RGBA &rgba);
	inline void setRed(real red)		{ c[0] = red;	}
	inline void setGreen(real green)	{ c[1] = green;	}
	inline void setBlue(real blue)		{ c[2] = blue;	}

	//////////////////////////////////////////////////////////////////////////
	// Index operators
	inline const real	operator[](unsigned int i) const { return c[i]; }
	inline real&		operator[](unsigned int i)		 { return c[i]; }

	inline const real r() const { return c[0];	}
	inline const real g() const { return c[1];	}
	inline const real b() const { return c[2];	}

	//////////////////////////////////////////////////////////////////////////
	// In place arithmetic
	inline Spectrum& operator+=(const Spectrum& rhs);
	inline Spectrum& operator*=(const Spectrum& rhs);
	inline Spectrum& operator/=(const Spectrum& rhs);
	inline Spectrum& operator*=(real rhs);
	inline Spectrum& operator/=(real rhs);

	inline const Spectrum operator+() const { return *this; }
	inline const Spectrum operator-() const;

	//////////////////////////////////////////////////////////////////////////
	// Manipulation
	const bool isBlack() const;
	void clamp(real min = 0, real max = 1);
	const Spectrum getClamped(real min = 0, real max = LARGE_i) const;
	const Spectrum& exponent();
	const Spectrum& sqrt();
	inline const real sum() const ;

	//////////////////////////////////////////////////////////////////////////
	// Static Methods
	static inline const Spectrum black()	{ return Spectrum(0.0f, 0.0f, 0.0f); }
	static inline const Spectrum red()		{ return Spectrum(1.0f, 0.0f, 0.0f); }
	static inline const Spectrum green()	{ return Spectrum(0.0f, 1.0f, 0.0f); }
	static inline const Spectrum blue() 	{ return Spectrum(0.0f, 0.0f, 1.0f); }
	static inline const Spectrum white()	{ return Spectrum(1.0f, 1.0f, 1.0f); }

	//////////////////////////////////////////////////////////////////////////
	// Friend Methods
	friend ostream& operator<<(ostream& out, const Spectrum& the_rgb);
	friend const Spectrum operator*(const Spectrum& c, real f);
	friend const Spectrum operator*(real f, const Spectrum& c);
	friend const Spectrum operator/(const Spectrum& c, real f);
	friend const Spectrum operator*(const Spectrum& c1, const Spectrum& c2);
	friend const Spectrum operator/(const Spectrum& c1, const Spectrum& c2);
	friend const Spectrum operator+(const Spectrum& c1, const Spectrum& c2);
	friend const Spectrum operator-(const Spectrum& c1, const Spectrum& c2);
	friend const Spectrum exp(const Spectrum&);
	friend const Spectrum log(const Spectrum&);

#ifdef _USE_XML_
	//////////////////////////////////////////////////////////////////////////
	// XML IO	
	static const std::string XmlElementName() { return std::string("Spectrum"); }
	void importXmlNode(xmlNodePtr node);
#endif

public :
	real c[N_COLOR_SAMPLES];
};

//////////////////////////////////////////////////////////////////////////
// Constructors
inline Spectrum::Spectrum(real v)						{ set(v);			}
inline Spectrum::Spectrum(real cs[N_COLOR_SAMPLES])		{ set(cs);			}
inline Spectrum::Spectrum(real c0, real c1, real c2)	{ set(c0, c1, c2);	}
inline Spectrum::Spectrum(const RGBA &rgba)				{ set(rgba);		}

//////////////////////////////////////////////////////////////////////////
// Assignment and set
inline void Spectrum::set(real v) {
	for(int i=0;i<N_COLOR_SAMPLES;++i)
		c[i] = v;
}

inline void Spectrum::set(real cs[N_COLOR_SAMPLES])	{ 
	std::memcpy(c, cs, N_COLOR_SAMPLES * sizeof(real));	
}

inline void Spectrum::set(real c0, real c1, real c2) { 
	c[0] = c0; c[1] = c1; c[2] = c2; 
}

inline void Spectrum::set(const Spectrum &rhs)	{ 
	(*this) = rhs;	
}

inline void Spectrum::set(const RGBA &rgba) {
	__ASSERT(N_COLOR_SAMPLES <= 4);
	for(int i=0;i<N_COLOR_SAMPLES;++i)
		c[i] = rgba[i];
}

//////////////////////////////////////////////////////////////////////////
// In place arithmetic
inline Spectrum& Spectrum::operator+=(const Spectrum& rhs)	{
	for(int i=0;i<N_COLOR_SAMPLES;++i)
		c[i] += rhs.c[i];
	return *this;
}

inline Spectrum& Spectrum::operator*=(const Spectrum& rhs)	{
	for(int i=0;i<N_COLOR_SAMPLES;++i)
		c[i] *= rhs.c[i];
	return *this;
}

inline Spectrum& Spectrum::operator/=(const Spectrum& rhs)	{
	for(int i=0;i<N_COLOR_SAMPLES;++i)
		c[i] /= rhs.c[i];
	return *this;
}

inline Spectrum& Spectrum::operator*=(real rhs)	{
	for(int i=0;i<N_COLOR_SAMPLES;++i)
		c[i] *= rhs;
	return *this;
}

inline Spectrum& Spectrum::operator/=(real rhs)	{
	for(int i=0;i<N_COLOR_SAMPLES;++i)
		c[i] /= rhs;
	return *this;
}

inline const Spectrum Spectrum::operator-() const { 
	Spectrum s;
	for(int i=0;i<N_COLOR_SAMPLES;++i)
		s.c[i] = -c[i];
	return s;
}	

//////////////////////////////////////////////////////////////////////////
// Manipulation
inline const real Spectrum::sum() const {
	real sum = 0;
	for(int i=0;i<N_COLOR_SAMPLES;++i)
		sum += c[i];
	return sum;
}

//////////////////////////////////////////////////////////////////////////
// Friend Methods
inline ostream& operator<<(ostream & os, const Spectrum& the_rgb)	{
	for(int i=0;i<N_COLOR_SAMPLES;++i)
		os << the_rgb.c[i] << ' ';
	return os;
}

inline const Spectrum operator*(const Spectrum& c, real f)	{	
	Spectrum s;
	for(int i=0;i<N_COLOR_SAMPLES;++i)
		s.c[i] = c.c[i] * f;
	return s;
}
inline const Spectrum operator*(real f, const Spectrum& c)	{	
	Spectrum s;
	for(int i=0;i<N_COLOR_SAMPLES;++i)
		s.c[i] = c.c[i] * f;
	return s;
}
inline const Spectrum operator/(const Spectrum& c, real f)	{	
	Spectrum s;
	for(int i=0;i<N_COLOR_SAMPLES;++i)
		s.c[i] = c.c[i] / f;
	return s;
}

inline const Spectrum operator*(const Spectrum& c1, const Spectrum& c2)	{	
	Spectrum s;
	for(int i=0;i<N_COLOR_SAMPLES;++i)
		s.c[i] = c1.c[i] * c2.c[i];
	return s;
}
inline const Spectrum operator/(const Spectrum& c1, const Spectrum& c2)	{	
	Spectrum s;
	for(int i=0;i<N_COLOR_SAMPLES;++i)
		s.c[i] = c1.c[i] / c2.c[i];
	return s;
}
inline const Spectrum operator+(const Spectrum& c1, const Spectrum& c2)	{	
	Spectrum s;
	for(int i=0;i<N_COLOR_SAMPLES;++i)
		s.c[i] = c1.c[i] + c2.c[i];
	return s;
}
inline const Spectrum operator-(const Spectrum& c1, const Spectrum& c2)	{	
	Spectrum s;
	for(int i=0;i<N_COLOR_SAMPLES;++i)
		s.c[i] = c1.c[i] - c2.c[i];
	return s;
}

inline const Spectrum exp(const Spectrum& c) { return Spectrum(::exp(c[0]), ::exp(c[1]), ::exp(c[2])); }
inline const Spectrum log(const Spectrum& c) { return Spectrum(::log(c[0]), ::log(c[1]), ::log(c[2])); }

#undef __ASSERT

};	// illuminationspace

#endif
