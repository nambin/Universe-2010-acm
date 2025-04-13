#include <stdafx.h>
#include <Nature/Fluid/CFDsolver/CFDIntegrator.h>
#include <Basis/Timer/Timer.h>

using namespace std;
using namespace basis;
using namespace geometryspace;
using namespace fluidspace;

CFDIntegrator_Liquid::CFDIntegrator_Liquid() : CFDIntegrator() {
	tracking = NULL;
	integrator = NULL;
	extrapolation = NULL;
	source = NULL;
	separation = NULL;
}

CFDIntegrator_Liquid::~CFDIntegrator_Liquid() {
	if(integrator) delete integrator;
	if(extrapolation) delete extrapolation;
	if(separation) delete separation;
}

void CFDIntegrator_Liquid::set(InterfaceTracking* _tracking, InterfaceIntegrator* _itg, Levelset_Extrapolation<real>* _ext, LiquidSource* s, GridBoundaryCondition_V_CeilingSeparation* sep) {
	tracking = _tracking;

	if(integrator)		delete integrator;		integrator = _itg;
	if(extrapolation)	delete extrapolation;	extrapolation = _ext;

	if(extrapolation)
		extrapolation->print_log(log_info);

	if(source)			delete source;			source = s;

	if(separation)		delete separation;		separation = sep;
	if(separation)		
		separation->print_log(log_info);
}

void CFDIntegrator_Liquid::update_one_time_step(CFDsolver& cfd, real dt) {
	if(cfd.getDimension() == 2) update_one_time_step_template((CFDsolver2D&)cfd, dt);
	else						update_one_time_step_template((CFDsolver3D&)cfd, dt);
}

void CFDIntegrator_Liquid::extrapolateVelocity(CFDsolver2D& cfd) {
	if(extrapolation == NULL)
		return;

	Grid2r phi2D; 
	tracking->getPhi(phi2D);

	extrapolation->extrapolate(cfd.u, MacGrid2D<real>::TYPE_U, *velBdry->uBdry, phi2D);
	extrapolation->extrapolate(cfd.v, MacGrid2D<real>::TYPE_V, *velBdry->vBdry, phi2D);
}

void CFDIntegrator_Liquid::extrapolateVelocity(CFDsolver3D& cfd) {
	if(extrapolation == NULL)
		return;

	Grid3r phi3D; 
	tracking->getPhi(phi3D);

	extrapolation->extrapolate(cfd.u, MacGrid3D<real>::TYPE_U, *velBdry->uBdry, phi3D);
	extrapolation->extrapolate(cfd.v, MacGrid3D<real>::TYPE_V, *velBdry->vBdry, phi3D);
	extrapolation->extrapolate(cfd.w, MacGrid3D<real>::TYPE_W, *velBdry->wBdry, phi3D);
}

void CFDIntegrator_Liquid::update_InternalBoundary(const FluidScene_Objects* objs, const CFDsolver& cfd) {
	objects = objs;
	integrator->update_InternalBoundary(objs, *tracking, cfd);
}

void CFDIntegrator_Liquid::importBinary(std::ifstream& ifs) {
	CFDIntegrator::importBinary(ifs);
	integrator->importBinary(ifs);
}

void CFDIntegrator_Liquid::exportBinary(std::ofstream& ofs) const {
	CFDIntegrator::exportBinary(ofs);
	integrator->exportBinary(ofs);
}

void CFDIntegrator_Liquid::post_process() const {
	CFDIntegrator::post_process();
	if(integrator)		integrator->post_process();
	if(extrapolation)	extrapolation->post_process();
}