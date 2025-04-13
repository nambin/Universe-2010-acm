#ifndef _VIRTUALSPACE_SCENELIGHT_H_
#define _VIRTUALSPACE_SCENELIGHT_H_

#include "SceneStateNode.h"
#include "../Light/Light.h"

namespace virtualspace {

using virtualspace::Light;

class SceneLight;
class PRINCIPLEDLL SceneLight_IOHandler : public IOHandler {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	SceneLight_IOHandler(SceneLight *);
	~SceneLight_IOHandler();

	//////////////////////////////////////////////////////////////////////////
	// Handler
	void keyboardHandler(unsigned char key, int	x, int y);
	void mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y);
	void mouseMotionHandler(int x, int y);

private : 
	SceneLight *sceneLight;
};

class Scene;
class PRINCIPLEDLL SceneLight : public SceneStateNode {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	SceneLight();
	SceneLight(const Light &);
	virtual ~SceneLight();

	//////////////////////////////////////////////////////////////////////////
	// Setter
	void setLight(const Light &_light)	{ light = _light;	}
	void setVisible(const bool _v)		{ visible = _v;		}

	//////////////////////////////////////////////////////////////////////////
	// Getter
	const Light &getLight() const	{ return light;		}
	const bool isVisible() const	{ return visible;	}

	//////////////////////////////////////////////////////////////////////////
	// gl Operation
	void glBind() const;
	void glUnbind() const;
	void glDraw(float size = 0.1f) const;

public :
	//////////////////////////////////////////////////////////////////////////
	// Idle Function
	void idle();	

protected : 
	Light light;
	bool visible;
};

};	// namespace virtualspace

#endif