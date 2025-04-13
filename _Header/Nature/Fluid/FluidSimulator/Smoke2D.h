#ifndef _FLUIDSPACE_SMOKE_2D_H_
#define _FLUIDSPACE_SMOKE_2D_H_

#include "../CFDsolver/CFDsolver2D.h"
#include "../DensityTracking/DensityTracking.h"
#include "../CFDsolver/CFDIntegrator.h"
#include "../FluidSimulator/FluidSimulator.h"
#include "../../Nature.h"
#include "../../../Principle/VirtualSpace/IOHandler/IOHandler.h"
#include "../../../Basis/Timer/Timer.h"

namespace fluidspace {

using namespace arithmetic;

class Smoke2D;
class NATUREDLL Smoke2D_IOHandler : public IOHandler {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Smoke2D_IOHandler(Smoke2D *_sim);
	~Smoke2D_IOHandler();

	//////////////////////////////////////////////////////////////////////////
	// Handler
	void keyboardHandler(unsigned char key, int	x, int y);
	void mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y);
	void mouseMotionHandler(int x, int y);

private : 
	Smoke2D *sim;
};

class NATUREDLL Smoke2D : public FluidSimulator2D {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Smoke2D();
	Smoke2D(grid_index_type nX, grid_index_type nY, real h, real dt, real maxCFL);	
	~Smoke2D();

	void post_process() const;

	void setComputingOption(const ParallelComputable& rhs);

	//////////////////////////////////////////////////////////////////////////
	// Initialization
	void set(grid_index_type nX, grid_index_type nY, real h, real dt, real maxCFL);	
	void setCoefficients(real g=-9.81f, real buoyancy=15.0f, real rho=1.0f, bool noSlip = false);

	void setDensity(const Grid2r &den);	

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	void glBind() const;
	void glDraw() const;	

	//////////////////////////////////////////////////////////////////////////
	// Object Coupling
	void update_InternalBoundary(const FluidScene_Objects*);

	//////////////////////////////////////////////////////////////////////////
	// Update
	void update(real dt);	

	//////////////////////////////////////////////////////////////////////////
	// Initial View Setting
	const Vector3f viewEye() const { return Vector3f(0.5f*solver.Lx, 0.5f*solver.Ly, MAX(solver.Lx, solver.Ly) * 1.25f); }
	const Vector3f viewRef() const { return Vector3f(0.5f*solver.Lx, 0.5f*solver.Ly, 0.0f);	}
	const Vector3f viewUp() const  { return Vector3f(0.0f, 1.0f, 0.0f);						}
	const float viewTheta() const  { return 45.0f;											}

public :
	CFDsolver2D				solver;
	DensityTracking2D		densityField;

	CFDIntegrator_Smoke*	integrator;
	
	//////////////////////////////////////////////////////////////////////////
	// Drawing Options
	bool dBoundingBox;
	bool dDensity;
	bool dPressure;
	bool dKinematicBoundary;
	bool dVelocity;

	real scaleVelocity;

	//////////////////////////////////////////////////////////////////////////
	// Friends 
	friend Smoke2D_IOHandler;
};

};

#endif