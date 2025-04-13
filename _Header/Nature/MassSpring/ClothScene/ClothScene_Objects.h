#ifndef _MASS_SPRING_SPACE_DYNAMIC_OBJECT_SCENE_H_
#define _MASS_SPRING_SPACE_DYNAMIC_OBJECT_SCENE_H_

#include "../../Nature/DynamicObject.h"
#include "../Environment/Environments.h"
#include "../../Nature.h"
#include <vector>

namespace springspace {

using nature::DynamicObject;
using std::vector;

class NATUREDLL ClothScene_Objects {
public : 
	ClothScene_Objects() {}

	void push_DynamicObject(DynamicObject<Vector3<real>>);

	struct QueryInfo {
		real			phi, mu;
		Vector3<real>	normal, vel;
	};

	//////////////////////////////////////////////////////////////////////////
	// Library for Collision Detection & Responce
	void query(QueryInfo&, const Vector3<real>&) const;

	void glDraw() const;

public : 
	vector<DynamicObject<Vector3<real>>>	objects;
};

};		// springspace

#endif