#ifndef _NATURE_FLUIDSPACE_FLUIDSCENE_H_
#define _NATURE_FLUIDSPACE_FLUIDSCENE_H_

#include "../../Nature.h"
#include "../../../Principle/VirtualSpace/Scene/SceneObject.h"
#include "../../../Principle/VirtualSpace/IOHandler/IOHandler.h"
#include "../Environment/Environments.h"
#include "FluidSceneIntegrator.h"
#include "FluidScene_Objects.h"

namespace fluidspace {

using namespace basis;
using namespace virtualspace;

class FluidScene;
class NATUREDLL FluidScene_IOHandler : public IOHandler {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	FluidScene_IOHandler(FluidScene *_sim);
	~FluidScene_IOHandler();

	//////////////////////////////////////////////////////////////////////////
	// Handler
	void keyboardHandler(unsigned char key, int	x, int y);
	void mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y);
	void mouseMotionHandler(int x, int y);

private : 
	FluidScene*		sim;
	SceneObject*	selectedObject;
};

class NATUREDLL FluidScene : public SceneObject {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	FluidScene();
	FluidScene(real _dt);
	~FluidScene();

	void release();

	virtual void post_process() const;

	//////////////////////////////////////////////////////////////////////////
	// Setter
	void set(real _dt);
	void setEndFrame(unsigned int _endframe) { endFrame = _endframe; }

	void enableSaveBinary(size_t _frame = 1.0e+5)	{ saveBinaries = true;	saving_frame = _frame;	}
	void disableSaveBinary()						{ saveBinaries = false;							}

	void enableSaveRenderObjectsOnly(size_t _frame = 1.0e+5)	{ saveRenderObjectsOnly = saveBinaries = true; saving_frame_renderObject = _frame;	}
	void disableSaveRenderObjectsOnly()							{ saveRenderObjectsOnly = false;													}

	void enableSaveImage()	{ saveImages = true;	}
	void disableSaveImage() { saveImages = false;	}	

	//////////////////////////////////////////////////////////////////////////
	// glBind
	virtual void glBind() const;
	void glDraw() const;

	void resizeScreen(unsigned int width, unsigned int height);

	//////////////////////////////////////////////////////////////////////////
	// Idle Function
	void idle();

	//////////////////////////////////////////////////////////////////////////
	// Save & Load
	void saveBinary() const;
	void saveRenderObject() const;
	void saveImage() const;

	void loadBinary(size_t frame);

	//////////////////////////////////////////////////////////////////////////
	// Simulator Initialization
	void initFluidSimulator(FluidSimulator*);
	void initFluidScene_Objects(FluidScene_Objects*);
	void initFluidSceneIntegrator(FluidSceneIntegrator*);

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

	const Vector3f viewEye_Side() const { return Vector3f(1.16162, 0.722146, 2.35805);		}
	const Vector3f viewRef_Side() const { return Vector3f(0.306223, 0.391154, 0.263992);	}
	const Vector3f viewUp_Side() const  { return Vector3f(-0.0775126, 0.989164, -0.124687); }
	const float viewTheta_Side() const  { return 45;										}

public : 
	unsigned int	frame, endFrame;
	real			dt, time;

	FluidSceneIntegrator*	integrator;

	FluidSimulator*			fluid;	
	FluidScene_Objects*		objects;

	//////////////////////////////////////////////////////////////////////////
	// Option	
	bool	isPause;
	bool	saveImages;
	bool	saveBinaries;				size_t saving_frame;
	bool	saveRenderObjectsOnly;		size_t saving_frame_renderObject;

	mutable bool isBeingSimulated;
};

};		// fluidspace

#endif