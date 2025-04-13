#include <stdafx.h>
#include <Nature/Fluid/InterfaceTracking/InterfaceIntegrator_Levelset.h>

using namespace std;
using namespace fluidspace;

InterfaceIntegrator_Levelset::InterfaceIntegrator_Levelset() {
	advect = NULL;
	redistancing = NULL;
}

InterfaceIntegrator_Levelset::~InterfaceIntegrator_Levelset() {
	if(advect) delete advect;
	if(redistancing) delete redistancing;	
}

void InterfaceIntegrator_Levelset::set_Advection(CFDAdvection* func) {
	if(advect) delete advect;
	advect = func;

	advect->print_log(log_info);
}

void InterfaceIntegrator_Levelset::set_Redistancing(Levelset_Redistancing<real>* func) {
	if(redistancing) delete redistancing;
	redistancing = func;

	redistancing->print_log(log_info);
}

void InterfaceIntegrator_Levelset::update(InterfaceTracking& lv, CFDsolver& cfd, real dt) {
	if(lv.getDimension() != cfd.getDimension()) {
		cout << "Error in InterfaceIntegrator_Levelset::update(InterfaceTracking& lv, CFDsolver& cfd, real _dt)" << endl;
		exit(0);
	}

	if(lv.getDimension() == 2 && cfd.getDimension() == 2) {	update_lv((InterfaceTracking_Levelset2D&)lv, (CFDsolver2D&)cfd, dt); }
	if(lv.getDimension() == 3 && cfd.getDimension() == 3) {	update_lv((InterfaceTracking_Levelset3D&)lv, (CFDsolver3D&)cfd, dt); }

	++numUpdates;
}

void InterfaceIntegrator_Levelset::update(InterfaceTracking& lv, VelocityField<real>& vel, real dt) {
	if(lv.getDimension() != vel.getDimension()) {
		cout << "Error in InterfaceIntegrator_Levelset::update(InterfaceTracking& lv, VelocityField<real>& vel, real dt)" << endl;
		exit(0);
	}

	if(lv.getDimension() == 2 && vel.getDimension() == 2) {	update_lv((InterfaceTracking_Levelset2D&)lv, (VelocityField2D<real>&)vel, dt); }
	if(lv.getDimension() == 3 && vel.getDimension() == 3) {	update_lv((InterfaceTracking_Levelset3D&)lv, (VelocityField3D<real>&)vel, dt); }

	++numUpdates;
}

void InterfaceIntegrator_Levelset::update_InternalBoundary(const FluidScene_Objects* objs, const InterfaceTracking& lv, const CFDsolver& cfd) {
	if(internalBdry)
		delete internalBdry; 
	internalBdry = NULL;

	if(objs->getDimension() != lv.getDimension() || objs->getDimension() != cfd.getDimension())
		throwError("Error in InterfaceIntegrator_Levelset::update_InternalBoundary(const FluidScene_Objects* objs, const InterfaceTracking& lv, const CFDsolver& cfd)");

	if(objs->getDimension() == 2)	internalBdry = new InternalBoundary_Implicit2D(10);
	else							internalBdry = new InternalBoundary_Implicit3D(10);

	internalBdry->update_InternalBoundary(objs, lv, cfd);
}