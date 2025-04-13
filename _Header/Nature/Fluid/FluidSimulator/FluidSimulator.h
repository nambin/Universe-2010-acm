#ifndef _FLUIDSPACE_FLUID_SIMULATOR_H_
#define _FLUIDSPACE_FLUID_SIMULATOR_H_

#include "../../../Principle/VirtualSpace/Scene/SceneObject.h"
#include "../../../Basis/ParallelProgramming/ParallelComputable.h"
#include "../../../Basis/TBBProgramming/TBBProgramming.h"
#include "../../../Basis/Basis/Log.h"
#include "../../Nature.h"
#include "../Environment/Environments.h"
#include "../FluidScene/FluidScene_Objects.h"

namespace fluidspace {

using namespace basis;
using namespace virtualspace;

class FluidScene_Objects;

class NATUREDLL FluidSimulator : public SceneObject, public ParallelComputable {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	FluidSimulator();
	FluidSimulator(real _dt, real _cfl);
	virtual ~FluidSimulator() = 0;

	virtual void post_process() const { log_info << "Fluid Simulator _ Have no Post Process" << std::endl; }

	//////////////////////////////////////////////////////////////////////////
	// Setter
	void set(real _dt, real _cfl);
	void setEndFrame(unsigned int _endframe) { endFrame = _endframe; }

	void enableSaveBinary()	{ saveBinaries = true;		}
	void disableSaveBinary() { saveBinaries = false;	}
	
	void enableSaveImage()	{ saveImages = true;	}
	void disableSaveImage() { saveImages = false;	}

	//////////////////////////////////////////////////////////////////////////
	// glBind
	virtual void glBind() const = 0;

	//////////////////////////////////////////////////////////////////////////
	// Idle Function
	void idle();

	//////////////////////////////////////////////////////////////////////////
	// Save 
	virtual void loadBinary(size_t frame);
	virtual void saveBinary(size_t frame) const;
	virtual void saveRenderObject(size_t frame) const;

	void importBinary(std::ifstream& ifs);
	void exportBinary(std::ofstream& ofs) const;

	void saveImage() const;

	//////////////////////////////////////////////////////////////////////////
	// Object Coupling
	virtual void update_InternalBoundary(const FluidScene_Objects*) { basis::throwError("Not Implemented _ FluidSimulator::update_InternalBoundary()"); }

	//////////////////////////////////////////////////////////////////////////
	// Update	
	virtual void update(real dt) = 0;
	void updateOneStep();

	virtual const unsigned int getDimension() const = 0;

public : 
	string			image_filename;

	unsigned int	frame, endFrame;
	real			dt, dt0, time, maxCFL;

	const FluidScene_Objects* objects;

	//////////////////////////////////////////////////////////////////////////
	// Option	
	bool	isPause;
	bool	saveBinaries;	
	bool	saveImages;	
};

class NATUREDLL FluidSimulator2D : public FluidSimulator {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	FluidSimulator2D() : FluidSimulator() {}
	FluidSimulator2D(real _dt, real _cfl) : FluidSimulator(_dt, _cfl) {}
	virtual ~FluidSimulator2D() {}

	const unsigned int getDimension() const { return 2; }
};

class NATUREDLL FluidSimulator3D : public FluidSimulator {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	FluidSimulator3D() : FluidSimulator() {}
	FluidSimulator3D(real _dt, real _cfl) : FluidSimulator(_dt, _cfl) {}
	virtual ~FluidSimulator3D() {}

	const unsigned int getDimension() const { return 3; }
};

};		// namespace fluidspace

#endif
