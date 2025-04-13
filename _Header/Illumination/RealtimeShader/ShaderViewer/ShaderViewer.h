#ifndef _SHADERSPACE_SHADER_VIEWER_H_
#define _SHADERSPACE_SHADER_VIEWER_H_

#include "../../IlluminationSpace.h"
#include "../../../Principle/Virtualspace/IOHandler/IOHandler.h"
#include "../../../Principle/Virtualspace/Scene/SceneObject.h"
#include "../../../Principle/VirtualSpace/Scene/Scene.h"
#include "../cgShader/cgShader_VolumeRendering.h"

namespace shaderspace {

using namespace virtualspace;

class ShaderViewer;
class ILLUMINATIONDLL ShaderViewer_IOHandler : public IOHandler {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	ShaderViewer_IOHandler(ShaderViewer *_sim);
	~ShaderViewer_IOHandler();

	//////////////////////////////////////////////////////////////////////////
	// Handler
	void keyboardHandler(unsigned char key, int	x, int y);
	void mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y);
	void mouseMotionHandler(int x, int y);

private : 
	ShaderViewer *viewer;
};

class ILLUMINATIONDLL ShaderViewer : public SceneObject {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	ShaderViewer();
	~ShaderViewer();

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	void glBind() const;
	void glDraw() const;	

	void idle() {}
	void resizeScreen(unsigned int width, unsigned int height);

	//////////////////////////////////////////////////////////////////////////
	// Initial View Setting
	const Vector3f viewEye() const { return Vector3f(0.5f*vr.Lx, 0.5f*vr.Ly, MAX(vr.Lx, vr.Ly) * 1.42f); }
	const Vector3f viewRef() const { return Vector3f(0.5f*vr.Lx, 0.5f*vr.Ly, 0.0f);				}
	const Vector3f viewUp() const  { return Vector3f(0.0f, 1.0f, 0.0f);							}
	const float viewTheta() const  { return 45.0f;												}

public :
	cgShader_VolumeRendering vr;
};

};

#endif