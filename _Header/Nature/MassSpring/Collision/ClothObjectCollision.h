#ifndef _MASS_SPRING_SPACE_CLOTH_OBJECT_COLLISION_H_
#define _MASS_SPRING_SPACE_CLOTH_OBJECT_COLLISION_H_

#include "../../Nature.h"
#include "../Environment/Environments.h"
#include "../Core/QuadDynamicMesh.h"
#include "../Core/TriDynamicMesh.h"
#include "../ClothScene/ClothScene_Objects.h"

namespace springspace {

class NATUREDLL ClothObjectCollision {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	ClothObjectCollision() {}
	virtual ~ClothObjectCollision() {}

	//////////////////////////////////////////////////////////////////////////
	// Update	
	virtual void treatCollision(const ClothScene_Objects&, real dt, QuadDynamicMesh<real>&) = 0;
	virtual void treatCollision(const ClothScene_Objects&, real dt, TriDynamicMesh<real>&) = 0;	
};

class NATUREDLL Bridson03_Collision : public ClothObjectCollision {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Bridson03_Collision(real relax);	
	~Bridson03_Collision();

	//////////////////////////////////////////////////////////////////////////
	// Update	
	void treatCollision(const ClothScene_Objects&, real dt, QuadDynamicMesh<real>&);
	void treatCollision(const ClothScene_Objects&, real dt, TriDynamicMesh<real>&);

	void treatCollision(const ClothScene_Objects&, real dt, const Vec3Array<real>& pos0, Vec3Array<real>& vel, Constraints<real>&) const;

public : 
	real relaxTime;
};

};			// springspace

#endif