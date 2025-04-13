#pragma once

#include "../Environment/Environments.h"
#include "../../Nature.h"
#include "../../Fluid/FluidSimulator/FluidSimulator.h"
#include "../../../Basis/Basis/Log.h"

namespace fluidspace {

using namespace basis;
using namespace virtualspace;

using arithmetic::Vector3f;

class InterfaceTracking_Viewer;
class NATUREDLL InterfaceTracking_Viewer_IOHandler : public IOHandler {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	InterfaceTracking_Viewer_IOHandler(InterfaceTracking_Viewer *_viewer);
	~InterfaceTracking_Viewer_IOHandler();

	//////////////////////////////////////////////////////////////////////////
	// Handler
	void keyboardHandler(unsigned char key, int	x, int y);
	void mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y);
	void mouseMotionHandler(int x, int y);

private : 
	InterfaceTracking_Viewer *viewer;
};

class NATUREDLL InterfaceTracking_Viewer : public FluidSimulator {
public : 
	InterfaceTracking_Viewer();	

	const bool load(const string& filename);
	const bool load_frame(int frame);
	const bool move_frame(int icr);

	void exportOBJ(const string&) const;

	void glBind() const;
	void glDraw() const;	

	void update(real dt);
	const unsigned int getDimension() const { return -1; }

	const bool load_startframe();
	const bool load_endframe();

	const Vector3f viewEye() const;
	const Vector3f viewRef() const;
	const Vector3f viewUp() const;
	const float viewTheta() const;

public : 
	InterfaceTracking*	tracking;

	std::string			filename;

	//////////////////////////////////////////////////////////////////////////
	// Drawing Options
	bool dBoundingBox;
	bool dLiquid;
	bool dGrid;	
	bool dPhi;
	bool dAuxiliary;
	bool dMesh;
	bool dInnerBdry;

	//////////////////////////////////////////////////////////////////////////
	// Friends 
	friend InterfaceTracking_Viewer_IOHandler;
};

};		// fluidspace