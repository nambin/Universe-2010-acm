#ifndef _ILLUMINATION_FILTER_H_
#define _ILLUMINATION_FILTER_H_

#include <string>
#include "../../../Basis/Basis/BasicLibrary.h"
#include "../../../Basis/XmlLibrary/XmlLibrary.h"
#include "../Illumination/Environments.h"
#include "../Illumination/Constants.h"

using namespace basis;

namespace illuminationspace {

class Filter {
public :
	Filter() {}	
	Filter(real xw, real yw) : xWidth(xw), yWidth(yw), invXWidth(1.0f/xw), invYWidth(1.0f/yw) {}
	
	virtual const real evaluate(real fx, real fy) const = 0;

#ifdef _USE_XML_
	//////////////////////////////////////////////////////////////////////////
	// XML IO
	static const std::string XmlElementName() { return std::string("Filter"); }
	virtual void importXmlNode(xmlNodePtr node);

	static Filter* createFromXmlNode(xmlNodePtr node);
#endif	

public :
	real xWidth, yWidth, invXWidth, invYWidth;
};

class BoxFilter : public Filter {
public :
	BoxFilter() {}	
	BoxFilter(real xw, real yw) : Filter(xw, yw) {}
	
	const real evaluate(real fx, real fy) const {
		return 1.0f;
	}	

#ifdef _USE_XML_
	BoxFilter(xmlNodePtr node) { importXmlNode(node); }

	//////////////////////////////////////////////////////////////////////////
	// XML IO
	static const std::string XmlElementName() { return std::string("BoxFilter"); }
	virtual void importXmlNode(xmlNodePtr node);
#endif
};

class TriangleFilter : public Filter {
public :
	TriangleFilter() {}	
	TriangleFilter(real xw, real yw) : Filter(xw, yw) {}
	
	const real evaluate(real fx, real fy) const {
		return MAX(static_cast<real>(0.0f), xWidth - abs(fx)) * MAX(static_cast<real>(0.0f), yWidth - abs(fy));
	}

#ifdef _USE_XML_
	TriangleFilter(xmlNodePtr node) { importXmlNode(node); }

	//////////////////////////////////////////////////////////////////////////
	// XML IO
	static const std::string XmlElementName() { return std::string("TriangleFilter"); }
	virtual void importXmlNode(xmlNodePtr node);
#endif
};

class GaussianFilter : public Filter {
public : 
	GaussianFilter() {}	
	GaussianFilter(real xw, real yw, real a) : Filter(xw, yw) {
		alpha = a;
		expX = ::exp(-alpha * xWidth * xWidth);
		expY = ::exp(-alpha * yWidth * yWidth);
	}

	const real evaluate(real x, real y) const {
		return Gaussian(x, expX) * Gaussian(y, expY);
	}
	const real Gaussian(real d, real expv) const {
		return MAX(static_cast<real>(0.0f), static_cast<real>(::exp(-alpha*d*d)-expv));
	}	

#ifdef _USE_XML_
	GaussianFilter(xmlNodePtr node) { importXmlNode(node); }

	//////////////////////////////////////////////////////////////////////////
	// XML IO
	static const std::string XmlElementName() { return std::string("GaussianFilter"); }
	virtual void importXmlNode(xmlNodePtr node);
#endif

public :
	real alpha, expX, expY;
};

class MitchellFilter : public Filter {
public:
	MitchellFilter() {}	
	MitchellFilter(real b, real c, real xw, real yw) : Filter(xw, yw) { 
		B = b; C = c; 
	}
	
	const real evaluate(real x, real y) const {
		return Mitchell1D(x * invXWidth) * Mitchell1D(y * invYWidth);
	}

	const real Mitchell1D(real x) const {
		x = abs(2.0f * x);
		if (x > 1.0f) return ((-B - 6*C) * x*x*x + (6*B + 30*C) * x*x + (-12*B - 48*C) * x + (8*B + 24*C)) * (1.f/6.f);
		else 		  return ((12 - 9*B - 6*C) * x*x*x + (-18 + 12*B + 6*C) * x*x + (6 - 2*B)) * (1.f/6.f);
	}	

#ifdef _USE_XML_
	MitchellFilter(xmlNodePtr node) { importXmlNode(node); }

	//////////////////////////////////////////////////////////////////////////
	// XML IO
	static const std::string XmlElementName() { return std::string("MitchellFilter"); }
	virtual void importXmlNode(xmlNodePtr node);
#endif

public :
	real B, C;
};

class LanczosSincFilter : public Filter {
public:
	LanczosSincFilter() {}	
	LanczosSincFilter(real xw, real yw,	real t) : Filter(xw, yw) {
		tau = t;
	}
	
	const real evaluate(real x, real y) const {
		return Sinc1D(x * invXWidth) * Sinc1D(y * invYWidth);
	}
	const real Sinc1D(real x) const {
		x = abs(x);
		if (x < 1e-5) return 1.f;
		if (x > 1.)   return 0.f;
		x *= PI_i;
		real sinc = sin(x * tau) / (x * tau);
		real lanczos = sin(x) / x;
		return sinc * lanczos;
	}	

#ifdef _USE_XML_
	LanczosSincFilter(xmlNodePtr node) { importXmlNode(node); }

	//////////////////////////////////////////////////////////////////////////
	// XML IO
	static const std::string XmlElementName() { return std::string("LanczosSincFilter"); }
	virtual void importXmlNode(xmlNodePtr node);
#endif

public :
	real tau;
};

};		// illuminationspace

#endif

