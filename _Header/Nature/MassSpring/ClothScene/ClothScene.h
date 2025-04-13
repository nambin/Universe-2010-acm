#ifndef _MASS_SPRING_SPACE_CLOTH_SCENE_H_
#define _MASS_SPRING_SPACE_CLOTH_SCENE_H_

#include "../../Nature.h"
#include "../../../Principle/VirtualSpace/Scene/SceneObject.h"
#include "../Environment/Environments.h"
#include "../Core/TriDynamicMesh.h"
#include "../Core/QuadDynamicMesh.h"
#include "../ClothScene/ClothScene_Objects.h"
#include "../ClothScene/ClothSceneIntegrator.h"

namespace springspace {

using namespace basis;
using namespace virtualspace;

class ClothScene;
class NATUREDLL ClothScene_IOHandler : public IOHandler {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	ClothScene_IOHandler(ClothScene *_sim);
	~ClothScene_IOHandler();

	//////////////////////////////////////////////////////////////////////////
	// Handler
	void keyboardHandler(unsigned char key, int	x, int y);
	void mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y);
	void mouseMotionHandler(int x, int y);

private : 
	ClothScene *sim;
};

class NATUREDLL ClothScene : public SceneObject {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	ClothScene();
	ClothScene(real _dt);
	~ClothScene();

	void release();

	//////////////////////////////////////////////////////////////////////////
	// Setter
	void set(real _dt);
	void setEndFrame(size_t _endframe) { endFrame = _endframe; }

	void enableSaveImage()	{ saveImages = true;	}
	void disableSaveImage() { saveImages = false;	}

	//////////////////////////////////////////////////////////////////////////
	// glBind
	virtual void glBind() const;
	void glDraw() const;

	//////////////////////////////////////////////////////////////////////////
	// Idle Function
	void idle();

	//////////////////////////////////////////////////////////////////////////
	// Save Images
	void saveImage() const;

	//////////////////////////////////////////////////////////////////////////
	// Simulator Initialization
	void initDynamicMesh(DynamicMesh*);
	void initClothScene_Objects(ClothScene_Objects*);
	void initClothSceneIntegrator(ClothSceneIntegrator*);

	//////////////////////////////////////////////////////////////////////////
	// Update	
	void update(real dt);
	void updateOneStep();

	//////////////////////////////////////////////////////////////////////////
	// Initial View Setting
	const Vector3f viewEye() const;
	const Vector3f viewRef() const;
	const Vector3f viewUp() const;
	const float viewTheta() const;

public : 
	size_t						frame, endFrame;
	real						dt, time;

	DynamicMesh*				mesh;	
	ClothScene_Objects*			objects;
	ClothSceneIntegrator*		integrator;

	//////////////////////////////////////////////////////////////////////////
	// Option	
	bool	isPause;
	bool	saveImages;

	//////////////////////////////////////////////////////////////////////////
	// Draw Option	
	bool	dSurface;
	bool	dParticle;
	bool	dStretchSpring;
	bool	dShearSpring;
	bool	dBendingSpring;
	bool	dObjects;
};

};		// namespace springspace

#endif
