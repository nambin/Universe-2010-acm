#include <stdafx.h>
#include <Nature/Fluid/InterfaceTracking/InterfaceIntegrator_Levelset.h>

using namespace std;
using namespace fluidspace;

InterfaceIntegrator_PLS::InterfaceIntegrator_PLS() : InterfaceIntegrator_Levelset() {	
	particleMnp = new PLS_ParticleMnp();
	particleAdvect = new PLS_PartlcleAdvector_RK2();
	errCorrector = new PLS_ErrorCorrection();
}

InterfaceIntegrator_PLS::~InterfaceIntegrator_PLS() {
	if(particleMnp) delete particleMnp;
	if(particleAdvect) delete particleAdvect;
	if(errCorrector) delete errCorrector;
}

void InterfaceIntegrator_PLS::initialize(InterfaceTracking& levelset) const {
	if(levelset.getDimension() == 2) {
		InterfaceTracking_PLS2D* pls = dynamic_cast<InterfaceTracking_PLS2D*>(&levelset);

		particleMnp->initialize(*pls);
		pls->printParticleInformation(log_info);
	}
	else if(levelset.getDimension() == 3) {
		InterfaceTracking_PLS3D* pls = dynamic_cast<InterfaceTracking_PLS3D*>(&levelset);

		particleMnp->initialize(*pls);
		pls->printParticleInformation(log_info);
	}	
}

void InterfaceIntegrator_PLS::update(InterfaceTracking& lv, CFDsolver& cfd, real dt) {
	if(lv.getDimension() != cfd.getDimension()) {
		cout << "Error in InterfaceIntegrator_PLS::update(InterfaceTracking& lv, CFDsolver& cfd, real dt)" << endl;
		exit(0);
	}

	if(lv.getDimension() == 2 && cfd.getDimension() == 2) { ++numUpdates; update_pls((InterfaceTracking_PLS2D&)lv, (CFDsolver2D&)cfd, dt); }
	if(lv.getDimension() == 3 && cfd.getDimension() == 3) { ++numUpdates; update_pls((InterfaceTracking_PLS3D&)lv, (CFDsolver3D&)cfd, dt); }
}


void InterfaceIntegrator_PLS::update(InterfaceTracking& lv, VelocityField<real>& vel, real dt) {
	if(lv.getDimension() != vel.getDimension()) {
		cout << "Error in InterfaceIntegrator_PLS::update(InterfaceTracking& levelset, VelocityField<real>& vel, real dt)" << endl;
		exit(0);
	}

	if(lv.getDimension() == 2 && vel.getDimension() == 2) {	++numUpdates; update_pls((InterfaceTracking_PLS2D&)lv, (VelocityField2D<real>&)vel, dt); }
	if(lv.getDimension() == 3 && vel.getDimension() == 3) {	++numUpdates; update_pls((InterfaceTracking_PLS3D&)lv, (VelocityField3D<real>&)vel, dt); }
}
