#pragma once

#include "../../../Principle/Geometry/DistanceFieldLibrary/DistanceFields.h"
#include "../../Nature.h"
#include "../Environment/Environments.h"

namespace fluidspace {

class NATUREDLL DensitySource {
public : 
	DensitySource() : time(0) {}

	virtual void add(DensityTracking*, real dt) const = 0;

	mutable real time;
};

class NATUREDLL DensitySource_DistanceField : public DensitySource {
public : 
	DensitySource_DistanceField(geometryspace::DistanceField<real>* f, real mag_per_sec, real t0 = 0, real t1 = 1.0e+6) 
	: field(f), magnitude(mag_per_sec), t_start(t0), t_end(t1) {}

	void add(DensityTracking*, real dt) const;

	geometryspace::DistanceField<real>*		field;
	real									magnitude, t_start, t_end;
};

};		// fluidspace