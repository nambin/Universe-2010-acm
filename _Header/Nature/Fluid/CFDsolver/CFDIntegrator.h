#ifndef _FLUIDSPACE_FLUID_INTEGRATOR_H_
#define _FLUIDSPACE_FLUID_INTEGRATOR_H_

#include "../../Nature.h"
#include "../FluidScene/FluidScene_Objects.h"
#include "../Environment/Environments.h"
#include "../CFDsolver/CFDsolver2D.h"
#include "../CFDsolver/CFDAddForce.h"
#include "../CFDsolver/CFDAdvection.h"
#include "../CFDsolver/CFDProjection.h"
#include "../CFDsolver/CFDVelocityBoundaryCondition.h"
#include "../CFDsolver/CFDTimeStepRestriction.h"
#include "../DensityTracking/DensitySource.h"
#include "../InterfaceTracking/LiquidSource.h"
#include "../InterfaceTracking/InterfaceIntegrator_Levelset.h"

namespace fluidspace {

class NATUREDLL CFDIntegrator : public ParallelComputable {
public : 
	CFDIntegrator(bool _use_variable_dt = true);
	virtual ~CFDIntegrator();

	void post_process() const;

	void setComputingOption(const ParallelComputable& rhs);

	void set_CFDAddForce(CFDAddForce*);
	void set_CFDAdvection(CFDAdvection*);
	void set_CFDProjection(CFDProjection*);
	void set_CFDVelocityBoundaryCondition(GridBdry*, GridBdry*, GridBdry*);
	void set_CFDTimeStepRestriction(CFDTimeStepRestriction*);

	void update(CFDsolver&, real dt);
	
	void update_substepping(CFDsolver&, real dt);
	void update_variable_timestep(CFDsolver&, real dt);
	
	virtual void update_one_time_step(CFDsolver&, real dt) = 0;
	virtual void update_InternalBoundary(const FluidScene_Objects*, const CFDsolver&) = 0;

	virtual void importBinary(std::ifstream& ifs);
	virtual void exportBinary(std::ofstream& ofs) const;

public : 
	const bool						use_variable_dt;

	CFDAddForce*					addforce;
	CFDAdvection*					advect;
	CFDProjection*					project;
	CFDVelocityBoundaryCondition*	velBdry;
	CFDTimeStepRestriction*			timeStepRestriction;

	const FluidScene_Objects*		objects;

	mutable size_t					numTimeSteps;
};

class NATUREDLL CFDIntegrator_Smoke : public CFDIntegrator {
public : 
	CFDIntegrator_Smoke();
	~CFDIntegrator_Smoke();

	void post_process() const;

	void set(DensityTracking*, CFDAdvection*, DensitySource* s = NULL);	

	void update_InternalBoundary(const FluidScene_Objects*, const CFDsolver&);

	void update_one_time_step(CFDsolver&, real dt);
	template<class CFD> void update_one_time_step_template(CFD&, real dt);

	void importBinary(std::ifstream& ifs);
	void exportBinary(std::ofstream& ofs) const;

public : 
	DensityTracking*		field;				// just pointing _ Not Delete
	CFDAdvection*			advect_density;
	DensitySource*			source;
};

class NATUREDLL CFDIntegrator_Liquid : public CFDIntegrator {
public : 
	CFDIntegrator_Liquid();
	~CFDIntegrator_Liquid();

	void post_process() const;

	void set(InterfaceTracking*, InterfaceIntegrator*, Levelset_Extrapolation<real>*, LiquidSource* s = NULL, GridBoundaryCondition_V_CeilingSeparation* sp = NULL);	

	void update_InternalBoundary(const FluidScene_Objects*, const CFDsolver&);

	void update_one_time_step(CFDsolver&, real dt);
	template<class CFD> void update_one_time_step_template(CFD&, real dt);

	void extrapolateVelocity(CFDsolver2D&);
	void extrapolateVelocity(CFDsolver3D&);	

	void importBinary(std::ifstream& ifs);
	void exportBinary(std::ofstream& ofs) const;
	
public : 
	InterfaceTracking*					tracking;				// just pointing _ Not Delete
	InterfaceIntegrator*				integrator;
	Levelset_Extrapolation<real>*		extrapolation;

	LiquidSource*						source;

	GridBoundaryCondition_V_CeilingSeparation* separation;
};

};		// fluidspace

#include "CFDIntegrator_Implementation.h"

#endif