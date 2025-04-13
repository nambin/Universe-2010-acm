#ifndef _FLUIDSPACE_LIQUID3D_H_
#define _FLUIDSPACE_LIQUID3D_H_

#include "../CFDsolver/CFDsolver3D.h"
#include "../InterfaceTracking/InterfaceTracking_PLS3D.h"
#include "../FluidSimulator/FluidSimulator.h"
#include "../CFDsolver/CFDIntegrator.h"
#include "../../Nature.h"
#include "../../../Principle/VirtualSpace/IOHandler/IOHandler.h"

namespace fluidspace {

using arithmetic::Vector3f;

class Liquid3D;
class NATUREDLL Liquid3D_IOHandler : public IOHandler {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Liquid3D_IOHandler(Liquid3D *_sim);
	~Liquid3D_IOHandler();

	//////////////////////////////////////////////////////////////////////////
	// Handler
	void keyboardHandler(unsigned char key, int	x, int y);
	void mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y);
	void mouseMotionHandler(int x, int y);

private : 
	Liquid3D *sim;
};

class NATUREDLL Liquid3D : public FluidSimulator3D {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Liquid3D();
	Liquid3D(grid_index_type nX, grid_index_type nY, grid_index_type nZ, real h, real dt, real maxCFL, real maxCFL_Interface);	
	~Liquid3D();

	void post_process() const;

	void release();

	//////////////////////////////////////////////////////////////////////////
	// Initialization
	void set(grid_index_type nX, grid_index_type nY, grid_index_type nZ, real h, real dt, real maxCFL, real maxCFL_Interface);	
	void setCoefficients(real _g=-9.81f, real _rhoPos=1.0f, real _rhoNeg=75.0f, bool noSlip = false, real narrowBand = 5.0f);	

	void setPhi(const Grid3r &_phi);
	void setPhi(const DistanceField3D<real>&);

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	void glBind() const;
	void glDraw() const;	

	void saveRenderObject(std::size_t frame) const;

	//////////////////////////////////////////////////////////////////////////
	// Object Coupling
	void update_InternalBoundary(const FluidScene_Objects*);

	//////////////////////////////////////////////////////////////////////////
	// Update
	void update(real _dt);	

	//////////////////////////////////////////////////////////////////////////
	// Initial View Setting
	const Vector3f viewEye() const { return tracking->viewEye();	}
	const Vector3f viewRef() const { return tracking->viewRef();	}
	const Vector3f viewUp() const  { return tracking->viewUp();		}
	const float viewTheta() const  { return tracking->viewTheta();	}

public :
	real							maxCFL_Interface;

	CFDsolver3D						solver;
	InterfaceTracking*				tracking;	

	CFDIntegrator_Liquid*			integrator;

	//////////////////////////////////////////////////////////////////////////
	// Drawing Options
	bool dBoundingBox;
	bool dSurface;
	bool dMesh;
	bool dInnerBoundary;

	//////////////////////////////////////////////////////////////////////////
	// Friends 
	friend Liquid3D_IOHandler;
};

};

#endif

