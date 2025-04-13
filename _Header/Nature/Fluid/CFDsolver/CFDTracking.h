#ifndef _FLUIDSPACE_CFD_TRACKING_H_
#define _FLUIDSPACE_CFD_TRACKING_H_

#include "../../Nature.h"
#include "../Environment/Environments.h"
#include "../CFDsolver/CFDsolver2D.h"
#include "../CFDsolver/CFDsolver3D.h"

namespace fluidspace {

using namespace arithmetic;

class NATUREDLL CFDTracking_EE {
public : 
	CFDTracking_EE() {}

	template<typename VEC, typename CFD>
	static const VEC forwardtrack(const VEC& pos, const CFD& cfd, real dt) { return pos + dt * cfd.getVel(pos); }

	template<typename VEC, typename CFD>
	static const VEC backtrack(const VEC& pos, const CFD& cfd, real dt) { return forwardtrack(pos, cfd, -dt); }

	template<typename VEC, typename VEL, typename BBOX>
	static const VEC forwardtrack(const VEC& pos, const VEL& vel, const BBOX& aabb, real dt) { return pos + dt * vel.getVel(pos); }

	template<typename VEC, typename VEL, typename BBOX>
	static const VEC backtrack(const VEC& pos, const VEL& vel, const BBOX& aabb, real dt) { return forwardtrack(pos, vel, aabb, -dt); }
};

class NATUREDLL CFDTracking_RK2 {
public : 
	CFDTracking_RK2() {}

	template<typename VEC, typename CFD>
	static const VEC forwardtrack(const VEC& pos, const CFD& cfd, real dt) { 		
		static const real half = 0.5;

		VEC mid = pos + (half * dt) * cfd.getVel(pos);
		cfd.constrainedBdry(mid);
		return pos + dt * cfd.getVel(mid);		
	}

	template<typename VEC, typename CFD>
	static const VEC backtrack(const VEC& pos, const CFD& cfd, real dt) { return forwardtrack(pos, cfd, -dt); }

	template<typename VEC, typename VEL, typename BBOX>
	static const VEC forwardtrack(const VEC& pos, const VEL& vel, const BBOX& aabb, real dt) { 
		static const real half = 0.5;

		VEC mid = pos + (half * dt) * vel.getVel(pos);
		aabb.constrainedBdry(mid);
		return pos + dt * vel.getVel(mid);	
	}

	template<typename VEC, typename VEL, typename BBOX>
	static const VEC backtrack(const VEC& pos, const VEL& vel, const BBOX& aabb, real dt) { return forwardtrack(pos, vel, aabb, -dt); }
};

class NATUREDLL CFDTracking_RK4 {
public : 
	CFDTracking_RK4() {}

	template<typename VEC, typename CFD>
	static const VEC forwardtrack(const VEC& pos, const CFD& cfd, real dt) { 	
		static const real half = 0.5;

		real half_dt = half * dt;
																	VEC v  = cfd.getVel(pos);
		VEC p1 = pos + half_dt * v;		cfd.constrainedBdry(p1);	VEC v1 = cfd.getVel(p1);
		VEC p2 = pos + half_dt * v1;	cfd.constrainedBdry(p2);	VEC v2 = cfd.getVel(p2);
		VEC p3 = pos + dt * v2;			cfd.constrainedBdry(p3);	VEC v3 = cfd.getVel(p3);

		return pos + (v + real(2) * v1 + real(2) * v2 + v3) * dt / real(6);
	}

	template<typename VEC, typename CFD>
	static const VEC backtrack(const VEC& pos, const CFD& cfd, real dt) { return forwardtrack(pos, cfd, -dt); }

	template<typename VEC, typename VEL, typename BBOX>
	static const VEC forwardtrack(const VEC& pos, const VEL& vel, const BBOX& aabb, real dt) { 
		static const real half = 0.5;

		real half_dt = half * dt;
																	VEC v  = vel.getVel(pos);
		VEC p1 = pos + half_dt * v;		aabb.constrainedBdry(p1);	VEC v1 = vel.getVel(p1);
		VEC p2 = pos + half_dt * v1;	aabb.constrainedBdry(p2);	VEC v2 = vel.getVel(p2);
		VEC p3 = pos + dt * v2;			aabb.constrainedBdry(p3);	VEC v3 = vel.getVel(p3);

		return pos + (v + real(2) * v1 + real(2) * v2 + v3) * dt / real(6);
	}

	template<typename VEC, typename VEL, typename BBOX>
	static const VEC backtrack(const VEC& pos, const VEL& vel, const BBOX& aabb, real dt) { return forwardtrack(pos, vel, aabb, -dt); }
};

};		// fluidspace

#endif
