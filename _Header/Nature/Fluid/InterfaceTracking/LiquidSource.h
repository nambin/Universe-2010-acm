#pragma once

#include "../../../Principle/Geometry/DistanceFieldLibrary/DistanceFields.h"
#include "../../Nature.h"
#include "../Environment/Environments.h"

namespace fluidspace {

class NATUREDLL LiquidSource {
public : 
	LiquidSource() : time(0) {}

	virtual void add(InterfaceTracking*, CFDsolver& cfd, real dt) const = 0;

//protected : 
	void setVel(CFDsolver2D&, const Vector2<real>&, const geometryspace::DistanceField2D<real>&, const real _eps = 1.5f) const;
	void setVel(CFDsolver3D&, const Vector3<real>&, const geometryspace::DistanceField3D<real>&, const real _eps = 1.5f) const;

	mutable real time;
};

class NATUREDLL LiquidSource_Continuous : public LiquidSource {
public : 
	LiquidSource_Continuous(geometryspace::DistanceField<real>* f) 
	: field(f) {}

	void add(InterfaceTracking*, CFDsolver& cfd, real dt) const;	
	
	void add_interval(real t0, real t1);

	void set_vel(const Vector2<real>& v) { vel2D = v; }
	void set_vel(const Vector3<real>& v) { vel3D = v; }

protected :
	Vector2<real>							vel2D;
	Vector3<real>							vel3D;

	geometryspace::DistanceField<real>*		field;
	std::vector<real>						interval;
};

class NATUREDLL LiquidSource_Discrete : public LiquidSource {
public : 
	LiquidSource_Discrete(geometryspace::DistanceField<real>* f) 
	: field(f) {}

	void add(InterfaceTracking*, CFDsolver& cfd, real dt) const;	

	void add_time_point(real t0);

	void set_vel(const Vector2<real>& v) { vel2D = v; }
	void set_vel(const Vector3<real>& v) { vel3D = v; }

protected : 
	Vector2<real>							vel2D;
	Vector3<real>							vel3D;

	geometryspace::DistanceField<real>*		field;	
	std::vector<real>						time_point;
};


};		// fluidspace