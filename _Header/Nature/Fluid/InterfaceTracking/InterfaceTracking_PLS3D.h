#ifndef _FLUIDSPACE_CFD_INTERFACE_TRACKING_PARTICLELEVELSET3D_H_
#define _FLUIDSPACE_CFD_INTERFACE_TRACKING_PARTICLELEVELSET3D_H_

#include "../../Nature.h"
#include <iostream>
#include "../Environment/Environments.h"
#include "../InterfaceTracking/InterfaceTracking_Levelset3D.h"

namespace fluidspace {

using geometryspace::MacGrid3D;
using std::vector;

class PLS_ParticleMnp;
class NATUREDLL InterfaceTracking_PLS3D : public InterfaceTracking_Levelset3D {
public : 
	typedef InterfaceTracking_Levelset3D::index_type index_type;

	static inline const real Epsilon() { return 1.0e-3; }

public : 
	InterfaceTracking_PLS3D() {}
	InterfaceTracking_PLS3D(index_type nX, index_type nY, index_type nZ, real h) : InterfaceTracking_Levelset3D(nX, nY, nZ, h) {}

	virtual void set(index_type nX, index_type nY, index_type nZ, real h);
	virtual void setPhi(const Grid3r&);
	virtual void setPhi(const geometryspace::DistanceField3D<real> &);

public :
	//////////////////////////////////////////////////////////////////////////
	// Utility
	template<typename STREAM> void printParticleInformation(STREAM &os = std::cout) const;

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	void glDrawAuxiliary() const { glDrawParticles(); }
	void glDrawParticles(const RGBA &color0 = RGBA::red(), const RGBA &color1 = RGBA::green(), GLfloat pointsize = 1.0f) const;

	const std::string name_class() const { return XmlElementName(); }
	static const std::string XmlElementName() { return std::string("InterfaceTracking_PLS3D"); }
	
public : 
	struct Particle {
		Vector3<real>	pos;
		char			sign;
		real			radius;
		bool			isEscaped;
		bool			isDeleted;

		Particle(const Vector3<real> &_pos, real _phi) 
		: pos(_pos), sign(SIGN(_phi)), radius(_phi * sign), isEscaped(false), isDeleted(false) {}

		const real getPhi(const Vector3<real> &_pos) const {
			return sign * (radius - length(pos,_pos));
		}
	};

public : 
	vector<Particle>	particles;
	vector<size_t>		idxEscapedParticles;	
};

template<typename STREAM> 
void InterfaceTracking_PLS3D::printParticleInformation(STREAM &os) const {

	const real eps = InterfaceTracking_PLS3D::Epsilon();

	size_t numPos = 0, numNeg = 0, numPosDeleted = 0, numNegDeleted = 0, numPosEscaped = 0, numNegEscaped = 0;
	size_t nearBoundary = 0;

	for(size_t i=0;i<particles.size();++i) {
		if(particles[i].sign == 1) {
			++numPos;

			if(particles[i].isDeleted) ++numPosDeleted;
			else if(particles[i].isEscaped) ++numPosEscaped;
		}
		else if(particles[i].sign == -1) {
			++numNeg;

			if(particles[i].isDeleted) ++numNegDeleted;
			else if(particles[i].isEscaped) ++numNegEscaped;
		}

		if(particles[i].pos[0] < eps || particles[i].pos[0] > Lx - eps 
			|| particles[i].pos[1] < eps || particles[i].pos[1] > Ly - eps
			|| particles[i].pos[2] < eps || particles[i].pos[2] > Lz - eps) ++nearBoundary;
	}

	os << "Number of Positive Particles : " << numPos << " _ Deleted : " << numPosDeleted << " _ Escaped : " << numPosEscaped << std::endl;
	os << "Number of Negative Particles : " << numNeg << " _ Deleted : " << numNegDeleted << " _ Escaped : " << numNegEscaped << std::endl;
	os << "Number of Particles near Boundary : " << nearBoundary << std::endl;
}

};		// fluidspace

#endif