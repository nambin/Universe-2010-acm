#include <stdafx.h>
#include <Nature/Fluid/CFDsolver/CFDAdvection.h>
#include <Principle/Geometry/MacGrid/MacGridFunctor.h>

using namespace std;
using namespace geometryspace;
using namespace fluidspace;

void CFDAdvection::advect(CFDsolver2D& cfd, real dt, const CFDVelocityBoundaryCondition& bdry) const {
	advect(cfd.u0, cfd.u, cfd, dt, MacGrid2D<real>::TYPE_U, *bdry.uBdry);
	advect(cfd.v0, cfd.v, cfd, dt, MacGrid2D<real>::TYPE_V, *bdry.vBdry);
	swap(cfd.u0, cfd.u);	swap(cfd.v0, cfd.v);
}

void CFDAdvection::advect(CFDsolver3D& cfd, real dt, const CFDVelocityBoundaryCondition& bdry) const {
	advect(cfd.u0, cfd.u, cfd, dt, MacGrid3D<real>::TYPE_U, *bdry.uBdry);
	advect(cfd.v0, cfd.v, cfd, dt, MacGrid3D<real>::TYPE_V, *bdry.vBdry);
	advect(cfd.w0, cfd.w, cfd, dt, MacGrid3D<real>::TYPE_W, *bdry.wBdry);
	swap(cfd.u, cfd.u0);	swap(cfd.v, cfd.v0);	swap(cfd.w, cfd.w0);	
}

void CFDAdvection::advect(DensityTracking2D& field, const CFDsolver2D& cfd, real dt) const {
	field.density.setNeumannBoundary();
	advect(field.density0, field.density, cfd, dt, MacGrid2D<real>::TYPE_S);
	swap(field.density0, field.density);
	field.density.setNeumannBoundary();
}

void CFDAdvection::advect(DensityTracking3D& field, const CFDsolver3D& cfd, real dt) const {
	field.getDensity().setNeumannBoundary();
	advect(field.getDensity0(), field.getDensity(), cfd, dt, MacGrid3D<real>::TYPE_S);
	swap(field.getDensity0(), field.getDensity());
	field.getDensity().setNeumannBoundary();
}

void CFDAdvection::advect(InterfaceTracking_Levelset2D& levelset, const CFDsolver2D& cfd, real dt) const {
	levelset.phi.setNeumannBoundary();
	advect(levelset.phi0, levelset.phi, cfd, dt, MacGrid2D<real>::TYPE_S);
	swap(levelset.phi0, levelset.phi);
	levelset.phi.setNeumannBoundary();	
}

void CFDAdvection::advect(InterfaceTracking_Levelset3D& levelset, const CFDsolver3D& cfd, real dt) const {
	levelset.phi.setNeumannBoundary();
	advect(levelset.phi0, levelset.phi, cfd, dt, MacGrid3D<real>::TYPE_S);
	swap(levelset.phi0, levelset.phi);
	levelset.phi.setNeumannBoundary();	
}

void CFDAdvection::advect(InterfaceTracking_Levelset2D& levelset, const VelocityField2D<real>& vel, real dt) const {
	levelset.phi.setNeumannBoundary();
	advect(levelset.phi0, levelset.phi, vel, dt);
	swap(levelset.phi0, levelset.phi);
	levelset.phi.setNeumannBoundary();	
}

void CFDAdvection::advect(InterfaceTracking_Levelset3D& levelset, const VelocityField3D<real>& vel, real dt) const {
	levelset.phi.setNeumannBoundary();
	advect(levelset.phi0, levelset.phi, vel, dt);
	swap(levelset.phi0, levelset.phi);
	levelset.phi.setNeumannBoundary();	
}