#ifndef _FLUIDSPACE_LIQUID2D_H_
#define _FLUIDSPACE_LIQUID2D_H_

#include "../../../Principle/VirtualSpace/IOHandler/IOHandler.h"
#include "../../Nature.h"
#include "../CFDsolver/CFDsolver2D.h"
#include "../InterfaceTracking/InterfaceTracking_PLS2D.h"
#include "../InterfaceTracking/InterfaceTracking_SRL2D.h"
#include "../InterfaceTracking/InterfaceIntegrator_SRL.h"
#include "../FluidSimulator/FluidSimulator.h"
#include "../CFDsolver/CFDIntegrator.h"

namespace fluidspace {

using arithmetic::Vector3f;

class Liquid2D;
class NATUREDLL Liquid2D_IOHandler : public IOHandler {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Liquid2D_IOHandler(Liquid2D *_sim);
	~Liquid2D_IOHandler();

	//////////////////////////////////////////////////////////////////////////
	// Handler
	void keyboardHandler(unsigned char key, int	x, int y);
	void mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y);
	void mouseMotionHandler(int x, int y);

private : 
	Liquid2D *sim;
};

class NATUREDLL Liquid2D : public FluidSimulator2D {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Liquid2D();
	Liquid2D(grid_index_type nX, grid_index_type nY, real h, real dt, real maxCFL, real maxCFL_Interface);
	~Liquid2D();

	void post_process() const;

	void release();

	//////////////////////////////////////////////////////////////////////////
	// Initialization
	void set(grid_index_type nX, grid_index_type nY, real h, real dt, real maxCFL, real maxCFL_Interface);
	void setCoefficients(real _g=-9.81f, real _rhoPos=1.0f, real _rhoNeg=75.0f, bool noSlip = false, real narrowBand = 5.0f);	

	void setPhi(const Grid2r &_phi);
	void setPhi(const DistanceField2D<real>&);

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	void glBind() const;
	void glDraw() const;	

	void loadBinary(std::size_t frame);
	void saveBinary(std::size_t frame) const;
	void saveRenderObject(std::size_t frame) const;

	//////////////////////////////////////////////////////////////////////////
	// Object Coupling
	void update_InternalBoundary(const FluidScene_Objects*);

	//////////////////////////////////////////////////////////////////////////
	// Update
	void update(real dt);

	//////////////////////////////////////////////////////////////////////////
	// Initial View Setting
	const Vector3f viewEye() const { return tracking->viewEye();	}
	const Vector3f viewRef() const { return tracking->viewRef();	}
	const Vector3f viewUp() const  { return tracking->viewUp();		}	
	const float viewTheta() const  { return tracking->viewTheta();	}

public :
	real							maxCFL_Interface;

	CFDsolver2D						solver;
	InterfaceTracking*				tracking;	

	CFDIntegrator_Liquid*			integrator;

	//////////////////////////////////////////////////////////////////////////
	// Drawing Options
	bool dBoundingBox;
	bool dInternalBoundary;
	bool dLiquid;
	bool dFilledInterface;
	bool dGrid;	
	bool dAuxiliary;
	bool dPhi;
	bool dParticle;
	bool dPressure;
	bool dVelocity;

	real scaleVelocity;

	//////////////////////////////////////////////////////////////////////////
	// Friends 
	friend Liquid2D_IOHandler;
};

};

#endif