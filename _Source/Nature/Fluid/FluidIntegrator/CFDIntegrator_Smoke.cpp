#include <stdafx.h>
#include <Nature/Fluid/CFDsolver/CFDIntegrator.h>
#include <Basis/Timer/Timer.h>

using namespace std;
using namespace basis;
using namespace geometryspace;
using namespace fluidspace;

CFDIntegrator_Smoke::CFDIntegrator_Smoke() : CFDIntegrator() { 
	field = NULL; 
	advect_density = NULL; 
	source = NULL;
}

CFDIntegrator_Smoke::~CFDIntegrator_Smoke() {
	if(advect_density) delete advect_density;
	if(source) delete source;
}

void CFDIntegrator_Smoke::set(DensityTracking* _field, CFDAdvection* _advect, DensitySource* _source) {
	if(_field)
		field = _field;

	if(_advect) {
		if(advect_density) delete advect_density;
		advect_density = _advect;

		log_info << "Density Advection _ "; _advect->print_log(log_info);
	}
	
	if(_source) {
		if(source) delete source;
		source = _source;
	}
}

void CFDIntegrator_Smoke::update_one_time_step(CFDsolver& cfd, real dt) {
	if(cfd.getDimension() == 2) update_one_time_step_template((CFDsolver2D&)cfd, dt);
	else						update_one_time_step_template((CFDsolver3D&)cfd, dt);
}

void CFDIntegrator_Smoke::update_InternalBoundary(const FluidScene_Objects* objs, const CFDsolver& cfd) {
	objects = objs;
}

void CFDIntegrator_Smoke::importBinary(std::ifstream& ifs) {
	CFDIntegrator::importBinary(ifs);
}

void CFDIntegrator_Smoke::exportBinary(std::ofstream& ofs) const {
	CFDIntegrator::exportBinary(ofs);
}

void CFDIntegrator_Smoke::post_process() const {
	CFDIntegrator::post_process();	
}