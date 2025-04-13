#include <stdafx.h>
#include <Nature/Fluid/InterfaceTracking/InterfaceIntegrator_SRL.h>

using namespace std;
using namespace fluidspace;

InterfaceIntegrator_SRL::~InterfaceIntegrator_SRL() {
	if(advect) delete advect;
	if(redistancing) delete redistancing;	
}

void InterfaceIntegrator_SRL::set_Advection(CFDAdvection* func) {
	if(advect) delete advect;
	advect = func;

	advect->print_log(log_info);
}

void InterfaceIntegrator_SRL::set_Redistancing(SRL_Redistancing* func) {
	if(redistancing) delete redistancing;
	redistancing = func;

	if(redistancing)	redistancing->print_log(log_info);
	else				log_info << "No re-distancing" << std::endl;
}

void InterfaceIntegrator_SRL::update(InterfaceTracking& srl, CFDsolver& cfd, real dt) {
	if(srl.getDimension() != cfd.getDimension()) {
		log_info << "Error in InterfaceIntegrator_SRL::update(InterfaceTracking& srl, VelocityField& vel, real _dt)" << endl;
		exit(0);
	}

	if(srl.getDimension() == 2 && cfd.getDimension() == 2) { update_srl((InterfaceTracking_SRL2D&)srl, (CFDsolver2D&)cfd, dt); }
	if(srl.getDimension() == 3 && cfd.getDimension() == 3) { update_srl((InterfaceTracking_SRL3D&)srl, (CFDsolver3D&)cfd, dt); }

	++numUpdates;
}

void InterfaceIntegrator_SRL::update(InterfaceTracking& srl, VelocityField<real>& vel, real dt) {
	if(srl.getDimension() != vel.getDimension()) {
		log_info << "Error in InterfaceIntegrator_SRL::update(InterfaceTracking& srl, VelocityField& vel, real _dt)" << endl;
		exit(0);
	}

	if(srl.getDimension() == 2 && vel.getDimension() == 2) { update_srl((InterfaceTracking_SRL2D&)srl, (VelocityField2D<real>&)vel, dt); }
	if(srl.getDimension() == 3 && vel.getDimension() == 3) { update_srl((InterfaceTracking_SRL3D&)srl, (VelocityField3D<real>&)vel, dt); }

	++numUpdates;
}

void InterfaceIntegrator_SRL::update_InternalBoundary(const FluidScene_Objects* objs, const InterfaceTracking& srl, const CFDsolver& cfd) {
	if(internalBdry)
		delete internalBdry; 
	internalBdry = NULL;

	if(objs->getDimension() != srl.getDimension() || objs->getDimension() != cfd.getDimension())
		throwError("Error in InterfaceIntegrator_SRL::update_InternalBoundary(const FluidScene_Objects* objs, const InterfaceTracking& srl, const CFDsolver& cfd)");

	if(objs->getDimension() == 2)	internalBdry = new InternalBoundary_Implicit2D(10);
	else							internalBdry = new InternalBoundary_Implicit3D(10);

	internalBdry->update_InternalBoundary(objs, srl, cfd);

	if(objs->getDimension() == 2)	((InterfaceTracking_SRL2D&)srl).update_InternalBoundary((InternalBoundary_Implicit2D*)internalBdry);
	else							((InterfaceTracking_SRL3D&)srl).update_InternalBoundary((InternalBoundary_Implicit3D*)internalBdry);
}

void InterfaceIntegrator_SRL::importBinary(std::ifstream& ifs) {
	if(redistancing)
		redistancing->importBinary(ifs);
}

void InterfaceIntegrator_SRL::exportBinary(std::ofstream& ofs) const {
	if(redistancing)
		redistancing->exportBinary(ofs);
}

void InterfaceIntegrator_SRL::post_process() const {
	if(advect) advect->post_process();
	if(redistancing) redistancing->post_process();

	log_info << "Time for SRL Tracking _ " << time_for_tracking << endl;
}