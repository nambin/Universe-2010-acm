#pragma once

#include "../../Nature.h"
#include "../Environment/Environments.h"
#include "../FluidScene/FluidScene_Objects.h"

namespace fluidspace {

class InterfaceTracking;
class InterfaceTracking_Levelset2D;
class InterfaceTracking_Levelset3D;
class InterfaceTracking_SRL2D;
class InterfaceTracking_SRL3D;
class CFDsolver;

class NATUREDLL InternalBoundary {
public : 
	InternalBoundary() : objects(NULL) {}

	virtual void update_InternalBoundary(const FluidScene_Objects* objs, const InterfaceTracking&, const CFDsolver&) { objects = objs; }	
	virtual void extrapolate_phi(InterfaceTracking&, real band) {
		throwError("Not Implemented _ InternalBoundary::update_InternalBoundary()");
	}

	virtual void glDraw() const { throwError("Not Implemented _ InternalBoundary::glDraw()"); }

	const FluidScene_Objects*	objects;
};

class NATUREDLL InternalBoundary_Implicit : public InternalBoundary {
public : 
	InternalBoundary_Implicit(real _band) : band(_band) {}

	void update_InternalBoundary(const FluidScene_Objects* objs, const InterfaceTracking&, const CFDsolver&) { 
		throwError("Not Implemented _ InternalBoundary_Implicit::update_InternalBoundary()");
	}	

	real band;
};

class NATUREDLL InternalBoundary_Implicit2D : public InternalBoundary_Implicit {
public : 
	InternalBoundary_Implicit2D(real _band) : InternalBoundary_Implicit(_band) {}

	void update_InternalBoundary(const FluidScene_Objects* objs, const InterfaceTracking&, const CFDsolver&);
	void extrapolate_phi(InterfaceTracking&, real band);

	void glDraw() const;

	Grid2r phi;

protected : 
	void extrapolate_phi(InterfaceTracking_Levelset2D&, real band);
	void extrapolate_phi(InterfaceTracking_SRL2D&, real band);
};

class NATUREDLL InternalBoundary_Implicit3D : public InternalBoundary_Implicit {
public : 
	InternalBoundary_Implicit3D(real _band) : InternalBoundary_Implicit(_band) {}

	void update_InternalBoundary(const FluidScene_Objects* objs, const InterfaceTracking&, const CFDsolver&);
	void extrapolate_phi(InterfaceTracking&, real band);

	void glDraw() const;

	Grid3r phi;

protected : 
	void extrapolate_phi(InterfaceTracking_Levelset3D&, real band);
	void extrapolate_phi(InterfaceTracking_SRL3D&, real band);
};

};		// fluidspace