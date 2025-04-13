#ifndef _FLUIDSPACE_PARTICLE_LEVELSET_PARTICLE_ADVECT_H_
#define _FLUIDSPACE_PARTICLE_LEVELSET_PARTICLE_ADVECT_H_

#include "../Environment/Random.h"
#include "../CFDsolver/CFDsolver2D.h"
#include "../CFDsolver/CFDsolver3D.h"
#include "../InterfaceTracking/InterfaceTracking_PLS2D.h"
#include "../InterfaceTracking/InterfaceTracking_PLS3D.h"

namespace fluidspace {

class NATUREDLL PLS_PartlcleAdvector {
public : 
	PLS_PartlcleAdvector() {}
	virtual void advect(InterfaceTracking_PLS2D&, const CFDsolver2D&, real dt) const = 0;
	virtual void advect(InterfaceTracking_PLS3D&, const CFDsolver3D&, real dt) const = 0;

	virtual void advect(InterfaceTracking_PLS2D&, const VelocityField2D<real>&, real dt) const = 0;
	virtual void advect(InterfaceTracking_PLS3D&, const VelocityField3D<real>&, real dt) const = 0;
};

class NATUREDLL PLS_PartlcleAdvector_RK2 : public PLS_PartlcleAdvector {
public : 
	PLS_PartlcleAdvector_RK2() {}
	void advect(InterfaceTracking_PLS2D&, const CFDsolver2D&, real dt) const;
	void advect(InterfaceTracking_PLS3D&, const CFDsolver3D&, real dt) const;

	void advect(InterfaceTracking_PLS2D&, const VelocityField2D<real>&, real dt) const;
	void advect(InterfaceTracking_PLS3D&, const VelocityField3D<real>&, real dt) const;
};

};		// fluidspace

#endif