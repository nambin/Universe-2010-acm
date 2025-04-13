#ifndef _ILLUMINATION_RAY_H_
#define _ILLUMINATION_RAY_H_

#include "../../../Principle/Arithmetic/Vector/Vector3.h"
#include "Spectrum.h"
#include "../../IlluminationSpace.h"
#include "../Illumination/Constants.h"

namespace illuminationspace {

using arithmetic::Vector3;

class ILLUMINATIONDLL Ray {
public :
	static inline const real Epsilon()	{ return 1.0e-3f;	}
	static inline const real Large()	{ return LARGE_i;	}

public :
	//////////////////////////////////////////////////////////////////////////
	// Constructors	
	inline Ray();
	inline Ray(const Vector3<real>& a, const Vector3<real>& b, const Spectrum& ex, const int depth);

	//////////////////////////////////////////////////////////////////////////
	// Setter
	inline void setRay(const Vector3<real>& a, const Vector3<real>& b);
	inline void setRay(const Vector3<real>& a, const Vector3<real>& b, const Spectrum& ex, const int depth);

	void setOrigin(const Vector3<real>& v)		{ data[0] = v;			}
	void setDirection(const Vector3<real>& v);	
	void setTmin(real _tmin) const				{ tmin = _tmin;			}
	void setTmax(real _tmax) const				{ tmax = _tmax;			}
	void setExtinction(const Spectrum& _ex)		{ _extinction = _ex;	}

	//////////////////////////////////////////////////////////////////////////
	// Getter
	Vector3<real>& origin()						{ return data[0];		}
	Vector3<real>& direction()					{ return data[1];		}

	const Vector3<real>& origin() const			{ return data[0];		}
	const Vector3<real>& direction() const		{ return data[1];		}
	const Vector3<real>& invDirection() const	{ return data[2];		}	
	const real tMin() const						{ return tmin;			}
	const real tMax() const						{ return tmax;			}
	const Spectrum& extinction() const			{ return _extinction;	}

	const Vector3<real> pointAtParameter(real t) const	{ return data[0] + t * data[1]; }
	const Vector3<real> operator()(real t) const		{ return data[0] + t * data[1]; }

public :
	Vector3<real>	data[3];
	int				posneg[3];

	mutable	int		depth;
	mutable real	tmin, tmax;
	Spectrum		_extinction;	

	// For Motion Blur
	real			time;
};

inline Ray::Ray() : depth(0), tmin(Ray::Epsilon()), tmax(Ray::Large()), _extinction(Spectrum(0,0,0)), time(0) {}

inline Ray::Ray(const Vector3<real>& a, const Vector3<real>& b, const Spectrum &_ex, const int _depth) 
: depth(0), tmin(Ray::Epsilon()), tmax(Ray::Large()), _extinction(Spectrum(0,0,0)), time(0) 
{ 
	setRay(a,b, _ex, _depth);	
}

inline void Ray::setRay(const Vector3<real>& a, const Vector3<real>& b)	{ 
	setRay(a,b,Spectrum(0.0f),0);	
}

inline void Ray::setRay(const Vector3<real>& a, const Vector3<real>& b, const Spectrum& _ex, const int _depth)	{ 
	data[0] = a; 
	setDirection(b);		
	
	depth = _depth;	
	tmin = Ray::Epsilon();	tmax = Ray::Large();
	_extinction = _ex; 

	time = 0;
}

inline std::ostream& operator<<(std::ostream& os, const Ray& r) {
	os << "(" << r.origin() << ") + t (" << r.direction() << ")";
	return os;
}

};		// illuminationspace

#endif