#ifndef _MASS_SPRING_SPACE_CLOTH_SELF_COLLISION_H_
#define _MASS_SPRING_SPACE_CLOTH_SELF_COLLISION_H_

#include "../../Nature.h"
#include "../../../Principle/Geometry/BVH/BVH.h"
#include "../Environment/Environments.h"
#include "../Core/QuadDynamicMesh.h"
#include "../Core/TriDynamicMesh.h"

namespace springspace {

using geometryspace::BVH;
using std::vector;
using std::pair;

class NATUREDLL ClothSelfCollision {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	ClothSelfCollision() {}
	virtual ~ClothSelfCollision() = 0 {}

	//////////////////////////////////////////////////////////////////////////
	// Update	
	virtual void treatCollision(real dt, QuadDynamicMesh<real>&) = 0;
	virtual void treatCollision(real dt, TriDynamicMesh<real>&) = 0;	
};

};			// springspace

#endif