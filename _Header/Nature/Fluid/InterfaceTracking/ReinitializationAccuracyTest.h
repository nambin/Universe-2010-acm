#pragma once

#include "../CFDsolver/CFDsolver2D.h"
#include "../InterfaceTracking/InterfaceTracking_PLS2D.h"
#include "../InterfaceTracking/InterfaceTracking_SRL2D.h"
#include "../FluidScene/FluidScene.h"
#include "../CFDsolver/CFDIntegrator.h"
#include "../../Nature.h"
#include "../../../Principle/VirtualSpace/IOHandler/IOHandler.h"
#include "SRL_Redistancing.h"

namespace fluidspace {

using arithmetic::Vector3f;

class ReinitializationAccuracyTest;
class NATUREDLL ReinitializationAccuracyTest_IOHandler : public IOHandler {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	ReinitializationAccuracyTest_IOHandler(ReinitializationAccuracyTest *_sim);
	~ReinitializationAccuracyTest_IOHandler();

	//////////////////////////////////////////////////////////////////////////
	// Handler
	void keyboardHandler(unsigned char key, int	x, int y);
	void mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y);
	void mouseMotionHandler(int x, int y);

private : 
	ReinitializationAccuracyTest *sim;
};

class NATUREDLL ReinitializationAccuracyTest : public FluidSimulator {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	ReinitializationAccuracyTest();
	~ReinitializationAccuracyTest() { release(); }

	void release();
	void update(real dt) {}

	//////////////////////////////////////////////////////////////////////////
	// Initialization
	void setExperiment_Russo();

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	void glBind() const;
	void glDraw() const;	

	const unsigned int getDimension() const { return -1; }

	//////////////////////////////////////////////////////////////////////////
	// Initial View Setting
	const Vector3f viewEye() const { if(tracking) return tracking->viewEye();	return Vector3f(0,0,1); }
	const Vector3f viewRef() const { if(tracking) return tracking->viewRef();	return Vector3f(0,0,0); }
	const Vector3f viewUp() const  { if(tracking) return tracking->viewUp();	return Vector3f(0,1,0); }
	const float viewTheta() const  { if(tracking) return tracking->viewTheta();	return 45.0f;			}


public :
	InterfaceTracking*				tracking;
	Levelset_Redistancing<real>*	redist;
	SRL_Redistancing*				redist_srl;

	//////////////////////////////////////////////////////////////////////////
	// Drawing Options
	bool dBoundingBox;
	bool dLiquid;
	bool dGrid;	
	bool dPhi;
	bool dAuxiliary;

	//////////////////////////////////////////////////////////////////////////
	// Friends 
	friend ReinitializationAccuracyTest_IOHandler;
};

};		// fluidspace

