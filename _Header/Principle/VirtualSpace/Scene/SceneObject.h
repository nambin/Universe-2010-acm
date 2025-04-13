#ifndef _VIRTUALSPACE_SCENEOBJECT_H_
#define _VIRTUALSPACE_SCENEOBJECT_H_

#include "ScenePrimitive.h"
#include "../Material/Material.h"
#include "../../Arithmetic/Vector/Vector3.h"

namespace virtualspace {

using namespace virtualspace;
using arithmetic::Vector3;

class Scene;
class PRINCIPLEDLL SceneObject : public ScenePrimitive {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	SceneObject(IOHandler *_io=NULL);
	virtual ~SceneObject();

	//////////////////////////////////////////////////////////////////////////
	// glBind
	virtual void glBind() const = 0;

	//////////////////////////////////////////////////////////////////////////
	// Idle Function
	virtual void idle() = 0;

	//////////////////////////////////////////////////////////////////////////
	// Material
	void setMaterial(const Material &mat) { material = mat; }
	const Material &getMaterial() const { return material;	}

	virtual const Vector3<float> viewEye() const = 0;
	virtual const Vector3<float> viewRef() const = 0;
	virtual const Vector3<float> viewUp() const = 0;
	virtual const float viewTheta() const = 0;

protected : 
	Material material;
};

};		// namespace virtualspace

#endif