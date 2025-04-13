#include <stdafx.h>
#include <Nature/Fluid/CFDsolver/CFDTimeStepRestriction.h>

using namespace std;
using namespace geometryspace;
using namespace fluidspace;

const real CFDTimeStepRestriction::getAllowedTimeStep(const CFDsolver& cfd, real h) const {
	const real C_cfl = cfd.getCFL(1.0f, h);

	log_info << "C_CFL : " << C_cfl << endl;

	if(IsAlmostZero(C_cfl))	return 1.0e+5;
	else					return maxCFL / C_cfl;
}

const size_t CFDTimeStepRestriction::getNumberOfSubStep(const CFDsolver& cfd, real dt, real h) const {
	const real C_cfl = cfd.getCFL(1.0f, h);
	const real cfl = dt * C_cfl;

	log_info << "CFL : " << cfl << endl;

	return static_cast<size_t>(cfl/maxCFL + 1);
}

const real CFDTimeStepRestriction::getAllowedTimeStep(const VelocityField<real>& vel, real h) const { 
	const real C_cfl = vel.getCFL(1.0f, h);

	if(IsAlmostZero(C_cfl))	return 1.0e+5;
	else					return maxCFL / C_cfl;
}

const size_t CFDTimeStepRestriction::getNumberOfSubStep(const VelocityField<real>& vel, real dt, real h) const { 
	const real C_cfl = vel.getCFL(1.0f, h);
	const real cfl = dt * C_cfl;

	return static_cast<size_t>(cfl/maxCFL + 1);
}

const real CFDTimeStepRestriction_SurfaceTension::getAllowedTimeStep(const CFDsolver& cfd, real h) const {
	const real C_cfl = cfd.getCFL(1.0f, h);
	const real G_cfl = sqrt(::abs(g) / h);
	const real S_cfl = sqrt(sigma / (MIN(rho_neg, rho_pos) * CB(h)));

	const real A_cfl = (C_cfl + sqrt( SQ(C_cfl) + 4 * SQ(G_cfl) + 4 * SQ(S_cfl) )) * 0.5f;	

	log_info << "A_CFL : " << A_cfl << endl;
	
	if(IsAlmostZero(A_cfl))	return 1.0e+5;
	else					return maxCFL / A_cfl;
}

const size_t CFDTimeStepRestriction_SurfaceTension::getNumberOfSubStep(const CFDsolver& cfd, real dt, real h) const {
	const real C_cfl = cfd.getCFL(1.0f, h);
	const real G_cfl = sqrt(::abs(g) / h);
	const real S_cfl = sqrt(sigma / (MIN(rho_neg, rho_pos) * CB(h)));

	const real A_cfl = (C_cfl + sqrt( SQ(C_cfl) + 4 * SQ(G_cfl) + 4 * SQ(S_cfl) )) * 0.5f;	
	const real cfl = A_cfl * dt;

	log_info << "CFL : " << cfl << endl;

	return static_cast<size_t>(cfl/maxCFL + 1);
}
