#ifndef _MASS_SPRING_DYNAMIC_MESH_H_
#define _MASS_SPRING_DYNAMIC_MESH_H_

#include "../../../Basis/Color/RGBA.h"
#include "../Core/ClothLibrary.h"
#include "../Core/Vec3Array.h"
#include "../Core/Constraints.h"
#include "../Core/MassSpringSystemMatrix.h"

namespace springspace {

class NATUREDLL DynamicMesh {
public : 
	typedef unsigned int index_type;

	DynamicMesh() {}

	virtual void postProcessing() = 0;

	virtual float width() const = 0;
	virtual float length() const = 0;
	virtual index_type size() const = 0;	

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	virtual void glDraw_FlatSurface(const Material& mf, const Material& mb) const {}
	virtual void glDraw_SmoothSurface(const Material& mf, const Material& mb) const {}

	virtual void glDraw_Particle(const RGBA& c0, const RGBA& cConstraint, GLfloat) const {}
	virtual void glDraw_StretchSpring(const RGBA&, GLfloat) const {}
	virtual void glDraw_ShearSpring(const RGBA&, GLfloat) const {}
	virtual void glDraw_BendingSpring(const RGBA&, GLfloat) const {}	
};

};			// springspace

#endif