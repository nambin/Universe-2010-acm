#pragma once

#include "../InterfaceTracking/InterfaceTracking_SLC3D.h"
#include "InterfaceIntegrator.h"

namespace fluidspace {

using namespace geometryspace;

class NATUREDLL InterfaceIntegrator_SLC : public InterfaceIntegrator {
public : 
	InterfaceIntegrator_SLC() {}
	void initialize(InterfaceTracking&) const {}
	
	void update(InterfaceTracking& slc, CFDsolver& cfd, real dt)			{ update_slc((InterfaceTracking_SLC3D&)slc, (CFDsolver3D&)cfd, dt);				}
	void update(InterfaceTracking& slc, VelocityField<real>& vel, real dt)	{ update_slc((InterfaceTracking_SLC3D&)slc, (VelocityField3D<real>&)vel, dt);	}

	template<class SLC, class CFD> void update_slc(SLC&, CFD&, real dt);

	virtual void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; }

	static const std::string XmlElementName() { return std::string("InterfaceIntegrator_SLC"); }
};

template<class CFD>
class VELField {
public : 
	VELField(const CFD& c) : cfd(c) {}
	const slcspace::SmVector3 operator()(const slcspace::SmVector3& _x) const {
		arithmetic::Vector3<double> vel = arithmetic::Vector3<double>(cfd.getVel(arithmetic::Vector3<real>(_x[0], _x[1], _x[2])));
		return slcspace::SmVector3(vel[0], vel[1], vel[2]);
	}
	const CFD&				cfd;
};

template<class VEL>
class TrackingSLC {
public : 
	TrackingSLC(slcspace::SlcSurface& _s, const VEL& v, const double _dt) : s(_s), vel(v), dt(_dt) {}
	const double operator()(const slcspace::SmVector3& _x) const {
		const slcspace::SmVector3 y = slcspace::traceBack(_x, vel, dt);
		return s.eval(y);
	}

	slcspace::SlcSurface&	s;
	const VEL&				vel;
	const double			dt;
};

template<class SLC, class CFD> 
void InterfaceIntegrator_SLC::update_slc(SLC& slc, CFD& cfd, real dt) {

	basis::Timer timer;

	TrackingSLC<VELField<CFD>> tracking_slc(*(slc.s), VELField<CFD>(cfd), dt);

	if(slc.s0)
		delete slc.s0;
	slc.s0 = new slcspace::SlcSurface();

	slc.s0->tree.buildTree(tracking_slc, slc.lc, slc.uc, slc.max_level);
	slc.s0->contourTree(tracking_slc);

	timer.logTime("Contouring SLC", log_info);	timer.reset();

	slc.s0->redistance();	

	timer.logTime("Redistacing SLC", log_info);	timer.reset();

	swap(slc.s0, slc.s);

	unsigned int numTriangles = slc.s->triangles.numTriangles();
	log_info << "Number of Triangles : " << numTriangles << endl;
}

};