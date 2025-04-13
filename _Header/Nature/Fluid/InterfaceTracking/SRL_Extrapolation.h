#pragma once

#include "../../Nature.h"
#include "InterfaceTracking_SRL2D.h"
#include "InterfaceTracking_SRL3D.h"

namespace fluidspace {

// SRL Extrapolation class for Internal Boundary
class NATUREDLL SRL_Extrapolation {
public : 
	typedef geometryspace::SRL<real>::index_type index_type;

	SRL_Extrapolation() {}

	virtual void extrapolate(InterfaceTracking_SRL2D&, const Grid2r&, const real zero_value) const { basis::throwError("Not Implemented _ SRL_Extrapolation::extrapolate(SRL2D<T>&)"); }	
	virtual void extrapolate(InterfaceTracking_SRL3D&, const Grid3r&, const real zero_value) const { basis::throwError("Not Implemented _ SRL_Extrapolation::extrapolate(SRL3D<T>&)"); }	
};

class NATUREDLL SRL_Extrapolation_PDE : public SRL_Extrapolation {
public : 
	SRL_Extrapolation_PDE() {}		

	void extrapolate(InterfaceTracking_SRL2D&, const Grid2r&, const real zero_value) const;
	void extrapolate(InterfaceTracking_SRL3D&, const Grid3r&, const real zero_value) const;
};

};			// fluidspace


