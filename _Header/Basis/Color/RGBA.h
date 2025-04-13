#ifndef _BASIS_RGBA_H_
#define _BASIS_RGBA_H_

#include <iostream>
#include <cassert>
#include "../Basis.h"
#include "../XmlLibrary/XmlLibrary.h"

namespace basis {

#ifdef _DEBUG
#define __ASSERT(c) assert(c)
#else
#define __ASSERT(c) {}
#endif

//////////////////////////////////////////////////////////////////////////
// RGBA Class Definition
//////////////////////////////////////////////////////////////////////////
class BASISDLL RGBA {
public : 
	typedef float value_type;

	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	inline RGBA();
	inline RGBA(const RGBA &rgba);
	inline RGBA(const RGBA::value_type c[4]);
	inline RGBA(RGBA::value_type r, RGBA::value_type g, RGBA::value_type b, RGBA::value_type a = 1.0f);

	//////////////////////////////////////////////////////////////////////////
	// Index operators
	inline const RGBA::value_type	operator[](std::size_t i) const;
	inline RGBA::value_type&		operator[](std::size_t i)		;

	//////////////////////////////////////////////////////////////////////////
	// Getter	
	inline const RGBA::value_type r() const { return data[0];	}
	inline const RGBA::value_type g() const { return data[1];	}
	inline const RGBA::value_type b() const { return data[2];	}
	inline const RGBA::value_type a() const { return data[3];	}

	//////////////////////////////////////////////////////////////////////////
	// Checker
	inline const bool isBlack() const;

	//////////////////////////////////////////////////////////////////////////
	// Get Pointer
	inline const RGBA::value_type* const getptr() const { return data; }

	//////////////////////////////////////////////////////////////////////////
	// Setter	
	inline void init();
	inline const RGBA& set(const RGBA &rgba);
	inline const RGBA& set(RGBA::value_type r, RGBA::value_type g, RGBA::value_type b, RGBA::value_type a);
	inline const RGBA& setRed(RGBA::value_type);
	inline const RGBA& setGreen(RGBA::value_type);
	inline const RGBA& setBlue(RGBA::value_type);
	inline const RGBA& setAlpha(RGBA::value_type);	

	inline const RGBA& operator=(const RGBA& rgba);

	//////////////////////////////////////////////////////////////////////////
	// In place arithmetic
	inline const RGBA& operator+=(const RGBA& rgba);
	inline const RGBA& operator*=(const RGBA& rgba);
	inline const RGBA& operator/=(const RGBA& rgba);
	inline const RGBA& operator*=(RGBA::value_type c);
	inline const RGBA& operator/=(RGBA::value_type c);

	inline const RGBA operator+() const;
	inline const RGBA operator-() const;

	//////////////////////////////////////////////////////////////////////////
	// Manipulation
	void clamp(RGBA::value_type min, RGBA::value_type max);

	//////////////////////////////////////////////////////////////////////////
	// Friend Methods
	inline friend std::ostream& operator<<(std::ostream& out, const RGBA& the_rgb);
	inline friend const RGBA operator*(const RGBA& c, RGBA::value_type f);
	inline friend const RGBA operator*(RGBA::value_type f, const RGBA& c);
	inline friend const RGBA operator/(const RGBA& c, RGBA::value_type f);
	inline friend const RGBA operator*(const RGBA& c1, const RGBA& c2);
	inline friend const RGBA operator/(const RGBA& c1, const RGBA& c2);
	inline friend const RGBA operator+(const RGBA& c1, const RGBA& c2);
	inline friend const RGBA operator-(const RGBA& c1, const RGBA& c2);

	//////////////////////////////////////////////////////////////////////////
	// Static Methods
	static const RGBA black()		{ return RGBA(0.0f, 0.0f, 0.0f, 1.0f); }
	static const RGBA red()			{ return RGBA(1.0f, 0.0f, 0.0f, 1.0f); }
	static const RGBA green()		{ return RGBA(0.0f, 1.0f, 0.0f, 1.0f); }
	static const RGBA blue() 		{ return RGBA(0.0f, 0.0f, 1.0f, 1.0f); }
	static const RGBA lightBlue()	{ return RGBA(173.0f/255.0f, 216.0f/255.0f, 230.0f/255.0f, 1.0f); }	
	static const RGBA pigmentBlue()	{ return RGBA(51.0f/255.0f, 51.0f/255.0f, 153.0f/255.0f, 1.0f); }	
	static const RGBA yellow() 		{ return RGBA(1.0f, 1.0f, 0.0f, 1.0f); }
	static const RGBA magenta()		{ return RGBA(1.0f, 0.0f, 1.0f, 1.0f); }	
	static const RGBA magentaDye()	{ return RGBA(202.0f/255.0f, 31.0f/255.0f, 123.0f/255.0f, 1.0f); }	
	static const RGBA brown()		{ return RGBA(150.0f/255.0f, 75.0f/255.0f, 0.0f, 1.0f); }	
	static const RGBA darkBrown()	{ return RGBA(101.0f/255.0f, 67.0f/255.0f, 33.0f/255.0f, 1.0f); }	
	static const RGBA paleBrown()	{ return RGBA(152.0f/255.0f, 118.0f/255.0f, 84.0f/255.0f, 1.0f); }	
	static const RGBA white()		{ return RGBA(1.0f, 1.0f, 1.0f, 1.0f); }

#ifdef _USE_XML_
	//////////////////////////////////////////////////////////////////////////
	// XML IO
	static inline const std::string XmlElementName() { return std::string("RGBA"); }
	void importXmlNode(xmlNodePtr node);
	const xmlNodePtr exportXmlNode() const;	
#endif

public : 
	RGBA::value_type data[4];
};

//////////////////////////////////////////////////////////////////////////
// Constructor/Destructor
inline RGBA::RGBA() {
	init();
}

inline RGBA::RGBA(const RGBA &rgba) {
	set(rgba);
}

inline RGBA::RGBA(const RGBA::value_type c[4]) {
	data[0] = c[0];	data[1] = c[1];	data[2] = c[2];	data[3] = c[3];
}

inline RGBA::RGBA(RGBA::value_type r, RGBA::value_type g, RGBA::value_type b, RGBA::value_type a) {
	data[0] = r; data[1] = g; data[2] = b; data[3] = a;
}

//////////////////////////////////////////////////////////////////////////
// Index operators
inline const RGBA::value_type RGBA::operator[](std::size_t i) const {
	__ASSERT(i < 4);
	return data[i];
}

inline RGBA::value_type& RGBA::operator[](std::size_t i)	{
	__ASSERT(i < 4);
	return data[i];
}

//////////////////////////////////////////////////////////////////////////
// Setter	
inline const bool RGBA::isBlack() const {
	if(data[0] == 0.0f && data[1] == 0.0f && data[2] == 0.0f)
		return true;
	else 
		return false;
}

//////////////////////////////////////////////////////////////////////////
// Setter	
inline void RGBA::init() {
	set(black());
}

inline const RGBA& RGBA::set(const RGBA &rgba) {
	for(std::size_t i=0;i<4;++i)
		data[i] = rgba[i];	
	return (*this);
}

inline const RGBA& RGBA::set(const RGBA::value_type r, const RGBA::value_type g, const RGBA::value_type b, const RGBA::value_type a) {
	data[0] = r; data[1] = g; data[2] = b; data[3] = a; return (*this);
}
inline const RGBA& RGBA::setRed(RGBA::value_type r)		{	data[0] = r;	return (*this); }
inline const RGBA& RGBA::setGreen(RGBA::value_type g)	{	data[1] = g;	return (*this); }
inline const RGBA& RGBA::setBlue(RGBA::value_type b)	{	data[2] = b;	return (*this); }
inline const RGBA& RGBA::setAlpha(RGBA::value_type a)	{	data[3] = a;	return (*this); }	

inline const RGBA& RGBA::operator=(const RGBA& rgba) {	return set(rgba);	}

//////////////////////////////////////////////////////////////////////////
// In place arithmetic
inline const RGBA& RGBA::operator+=(const RGBA& rgba) {
	for(std::size_t i=0;i<3;++i)
		data[i] += rgba[i];
	return (*this);
}

inline const RGBA& RGBA::operator*=(const RGBA& rgba) {
	for(std::size_t i=0;i<3;++i)
		data[i] *= rgba[i];
	return (*this);
}

inline const RGBA& RGBA::operator/=(const RGBA& rgba) {
	for(std::size_t i=0;i<3;++i)
		data[i] /= rgba[i];
	return (*this);
}

inline const RGBA& RGBA::operator*=(RGBA::value_type c) {
	for(std::size_t i=0;i<3;++i)
		data[i] *= c;
	return (*this);
}

inline const RGBA& RGBA::operator/=(RGBA::value_type c) {
	for(std::size_t i=0;i<3;++i)
		data[i] /= c;
	return (*this);
}

inline const RGBA RGBA::operator+() const {
	return (*this);
}

inline const RGBA RGBA::operator-() const {
	RGBA rgba = (*this) * (-1.0f);
	return rgba;
}

//////////////////////////////////////////////////////////////////////////
// Friend Methods
inline std::ostream& operator<<(std::ostream& os, const RGBA& rgba) {
	os << "[" << rgba[0] << ", " << rgba[1] << ", " << rgba[2] << ", " << rgba[3] << "]";
	return os;
}

inline const RGBA operator*(const RGBA& c, RGBA::value_type f) {
	RGBA rgba;
	for(int i=0;i<3;++i)
		rgba[i] = c[i] * f;
	rgba[3] = c[3];
	return rgba;
}

inline const RGBA operator*(RGBA::value_type f, const RGBA& c) {
	return c*f;
}

inline const RGBA operator/(const RGBA& c, RGBA::value_type f) {
	return c*(1.0f/f);
}

inline const RGBA operator*(const RGBA& c1, const RGBA& c2) {
	RGBA rgba;
	for(int i=0;i<3;++i)
		rgba[i] = c1[i] * c2[i];
	rgba[3] = c1[3];
	return rgba;
}

inline const RGBA operator/(const RGBA& c1, const RGBA& c2) {
	RGBA rgba;
	for(int i=0;i<3;++i)
		rgba[i] = c1[i] / c2[i];
	rgba[3] = c1[3];
	return rgba;
}

inline const RGBA operator+(const RGBA& c1, const RGBA& c2) {
	RGBA rgba;
	for(int i=0;i<3;++i)
		rgba[i] = c1[i] + c2[i];
	rgba[3] = c1[3];
	return rgba;
}

inline const RGBA operator-(const RGBA& c1, const RGBA& c2) {
	RGBA rgba;
	for(int i=0;i<3;++i)
		rgba[i] = c1[i] - c2[i];
	rgba[3] = c1[3];
	return rgba;
}

#undef __ASSERT

};		// namespace basis

#endif

