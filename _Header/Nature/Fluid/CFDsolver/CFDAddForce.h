#ifndef _FLUIDSPACE_CFD_ADD_FORCE_H_
#define _FLUIDSPACE_CFD_ADD_FORCE_H_

#include "../../Nature.h"
#include "../Environment/Environments.h"
#include "../CFDsolver/CFDsolver2D.h"
#include "../CFDsolver/CFDsolver3D.h"
#include "../DensityTracking/DensityTracking.h"
#include "../InterfaceTracking/InterfaceTracking.h"

namespace fluidspace {

class NATUREDLL CFDAddForce {
public : 
	typedef CFDsolver2D::index_type index_type;

	CFDAddForce() {}
	virtual void post_process() const { log_info << "CFDAddForce _ Have no Post Process" << std::endl; }

	virtual void addForce(CFDsolver2D&, real dt) = 0;
	virtual void addForce(CFDsolver3D&, real dt) = 0;

	void addForce(CFDsolver2D&, const Grid2<Vector2<real>,real>& field, real dt);
	void addForce(CFDsolver3D&, const Grid3<Vector3<real>,real>& field, real dt);

	virtual void print_log(Log& log_) const = 0;
	static const std::string XmlElementName() { return std::string("CFDAddForce"); }	
};

class NATUREDLL CFDAddForce_Gravity : public CFDAddForce {
public : 
	CFDAddForce_Gravity(real _g) : CFDAddForce(), g(_g) {}
	void addForce(CFDsolver2D&, real dt);
	void addForce(CFDsolver3D&, real dt);

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << " _ Gravity : " << g << endl; }
	static const std::string XmlElementName() { return std::string("CFDAddForce_Gravity"); }

	const real g;
};

class NATUREDLL CFDAddForce_DensityForce : public CFDAddForce_Gravity {
public : 
	CFDAddForce_DensityForce(real _g, real _buoyancy) : CFDAddForce_Gravity(_g), buoyancy(_buoyancy), density(NULL) {}
	void addForce(CFDsolver2D&, real dt);
	void addForce(CFDsolver3D&, real dt);

	void setDensityTracking(const DensityTracking* ptr) { density = ptr; }

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; CFDAddForce_Gravity::print_log(log_); }
	static const std::string XmlElementName() { return std::string("CFDAddForce_DensityForce"); }

	const DensityTracking*		density;		// just pointing _ Not Delete
	const real					buoyancy;	
};

class NATUREDLL CFDAddForce_SurfaceTension : public CFDAddForce_Gravity {
public : 
	CFDAddForce_SurfaceTension(real _g, real _sigma, real _rho_neg, real _rho_pos) 
	: CFDAddForce_Gravity(_g), sigma(_sigma), rho_neg(_rho_neg), rho_pos(_rho_pos), tracking(NULL) {}

	void set_InterfaceTracking(const InterfaceTracking* track) { tracking = track; }

	void addForce(CFDsolver2D&, real dt);
	void addForce(CFDsolver3D&, real dt);

	void print_log(Log& log_) const { 
		log_ << XmlElementName().c_str() << " _ Tension : " << sigma << " _ Density (air, water) : " << rho_pos << ", " << rho_neg << endl; 
		CFDAddForce_Gravity::print_log(log_); 
	}
	static const std::string XmlElementName() { return std::string("CFDAddForce_SurfaceTension"); }

	const InterfaceTracking*	tracking;		// just pointing _ Not Delete
	const real					sigma, rho_neg, rho_pos;
};

//// Not Using _ Weird
//class NATUREDLL CFDAddForce_Liquid_Bouyancy : public CFDAddForce_Gravity {
//public : 
//	CFDAddForce_Liquid_Bouyancy(real g, real rho_n, real rho_p) : CFDAddForce_Gravity(g), rho_neg(rho_n), rho_pos(rho_p) {}
//
//	void set_InterfaceTracking(const InterfaceTracking* track) { tracking = track; }
//
//	void addForce(CFDsolver2D&, real dt);
//	void addForce(CFDsolver3D&, real dt);
//
//	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; CFDAddForce_Gravity::print_log(log_); }
//	static const std::string XmlElementName() { return std::string("CFDAddForce_Liquid_Bouyancy"); }
//
//	const InterfaceTracking*	tracking;
//	real						rho_pos, rho_neg;
//};

};		// fluidspace

#endif