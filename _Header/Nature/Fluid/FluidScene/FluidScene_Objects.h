#ifndef _FLUID_SPACE_DYNAMIC_OBJECT_SCENE_H_
#define _FLUID_SPACE_DYNAMIC_OBJECT_SCENE_H_

#include "../../../Principle/VirtualSpace/Scene/SceneObject.h"
#include "../../Nature/DynamicObject.h"
#include "../../Nature.h"
#include <vector>

namespace fluidspace {

using namespace virtualspace;
using nature::DynamicObject;
using std::vector;

class FluidScene_Objects;
class NATUREDLL FluidScene_Objects_IOHandler : public IOHandler {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	FluidScene_Objects_IOHandler(FluidScene_Objects *_sim);
	~FluidScene_Objects_IOHandler();

	//////////////////////////////////////////////////////////////////////////
	// Handler
	void keyboardHandler(unsigned char key, int	x, int y);
	void mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y);
	void mouseMotionHandler(int x, int y);

private : 
	FluidScene_Objects*		objs;	
};

class NATUREDLL FluidScene_Objects : public SceneObject {
public : 
	FluidScene_Objects() : dObject(true) { ioHandler = new FluidScene_Objects_IOHandler(this); }
	virtual ~FluidScene_Objects() {}	

	virtual void post_process() const { log_info << "Fluid Scene Objects _ Have no Post Process" << std::endl; }

	void glBind() const { SceneObject::glBind(); glDraw(); } 
	void idle() {}

	virtual void glDraw() const = 0;

	virtual const unsigned int getDimension() const = 0;

	mutable bool dObject;

	const Vector3f viewEye() const { return Vector3f(0,0,1);	}
	const Vector3f viewRef() const { return Vector3f(0,0,0);	}
	const Vector3f viewUp() const  { return Vector3f(0,1,0);	}
	const float viewTheta() const  { return 45.0f;				}
};

class NATUREDLL FluidScene_Objects2D : public FluidScene_Objects {
public : 
	FluidScene_Objects2D() {}

	void push_DynamicObject(DynamicObject<Vector2<real>>);

	const bool isInside(const Vector2<real>&) const;

	const Vector2<real> getVel(const Vector2<real>& p) const;

	void glDraw() const;

	const unsigned int getDimension() const { return 2; }

public : 
	vector<DynamicObject<Vector2<real>>> objects;
};

class NATUREDLL FluidScene_Objects3D : public FluidScene_Objects {
public : 
	FluidScene_Objects3D() {}

	void push_DynamicObject(DynamicObject<Vector3<real>>);

	const bool isInside(const Vector3<real>&) const;

	const Vector3<real> getVel(const Vector3<real>& p) const;

	void glDraw() const;

	const unsigned int getDimension() const { return 3; }

public : 
	vector<DynamicObject<Vector3<real>>> objects;
};

};		// springspace

#endif