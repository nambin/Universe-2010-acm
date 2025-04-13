#ifndef _FLUIDSPACE_SMOKE_3D_H_
#define _FLUIDSPACE_SMOKE_3D_H_

#include "../../Nature.h"
#include "../CFDsolver/CFDsolver3D.h"
#include "../FluidSimulator/FluidSimulator.h"
#include "../DensityTracking/DensityTracking.h"
#include "../CFDsolver/CFDIntegrator.h"
#include "../../../Principle/VirtualSpace/IOHandler/IOHandler.h"
#include "../../../Basis/Timer/Timer.h"

namespace fluidspace {

using namespace arithmetic;

class Smoke3D;
class NATUREDLL Smoke3D_IOHandler : public IOHandler {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Smoke3D_IOHandler(Smoke3D *_sim);
	~Smoke3D_IOHandler();

	//////////////////////////////////////////////////////////////////////////
	// Handler
	void keyboardHandler(unsigned char key, int	x, int y);
	void mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y);
	void mouseMotionHandler(int x, int y);

private : 
	Smoke3D *sim;
};

class NATUREDLL Smoke3D : public FluidSimulator3D {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Smoke3D();
	Smoke3D(grid_index_type nX, grid_index_type nY, grid_index_type nZ, real h, real _dt, real _maxCFL);	
	~Smoke3D();

	void post_process() const;

	void setComputingOption(const ParallelComputable& rhs);

	//////////////////////////////////////////////////////////////////////////
	// Initialization
	void set(grid_index_type nX, grid_index_type nY, grid_index_type nZ, real h, real dt, real maxCFL);	
	void setCoefficients(real g=-9.81f, real buoyancy=13.0f, real rho=1.0f, bool noSlip = false);		
	
	void setDensity(const Grid3r &den);	

	//////////////////////////////////////////////////////////////////////////
	// gl Operation
	void glBind() const;
	void glDraw() const;	

	void resizeScreen(unsigned int width, unsigned int height);

	//////////////////////////////////////////////////////////////////////////
	// Object Coupling
	void update_InternalBoundary(const FluidScene_Objects*);

	//////////////////////////////////////////////////////////////////////////
	// Update
	void update(real _dt);

	//////////////////////////////////////////////////////////////////////////
	// Initial View Setting
	const Vector3f viewEye() const { return Vector3f(0.5f*solver.Lx, 0.5f*solver.Ly, MAX(solver.Lx, solver.Ly) * 2.25f); }
	const Vector3f viewRef() const { return Vector3f(0.5f*solver.Lx, 0.5f*solver.Ly, 0.0f);	}
	const Vector3f viewUp() const  { return Vector3f(0.0f, 1.0f, 0.0f);						}
	const float viewTheta() const  { return 45.0f;											}

public :
	CFDsolver3D			solver;
	DensityTracking3D	densityField;

	CFDIntegrator*		integrator;

	//////////////////////////////////////////////////////////////////////////
	// Variable for Parallel Programming	
	bool useAutoPartitioner;
	size_t grainSize[3];	

	//////////////////////////////////////////////////////////////////////////
	// Drawing Options
	bool dBoundingBox;
	bool dDensity;
	bool dVelocity;

	real scaleVelocity;

	//////////////////////////////////////////////////////////////////////////
	// Friends 
	friend Smoke3D_IOHandler;
};

};

#endif