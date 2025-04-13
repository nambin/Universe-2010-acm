#ifndef _MASS_SPRING_SPACE_PROVOT_LIMITING_H_
#define _MASS_SPRING_SPACE_PROVOT_LIMITING_H_

#include <vector>
#include "../../../Basis/Basis/Array.h"
#include "../../../Principle/Arithmetic/Vector/Vector3.h"
#include "../Environment/Environments.h"
#include "../Core/Spring.h"
#include "../Core/Constraints.h"
#include "../Core/TriDynamicMesh.h"
#include "../Core/QuadDynamicMesh.h"
#include "../../Nature.h"

namespace springspace {

using std::vector;
using basis::Array;
using arithmetic::Vector3;

class NATUREDLL PostProcessingDynamics {
public : 
	virtual void postProcess(real dt, TriDynamicMesh<real>&) const = 0;
	virtual void postProcess(real dt, QuadDynamicMesh<real>&) const = 0;
};

class NATUREDLL ProvotLimiting : public PostProcessingDynamics {
public : 
	ProvotLimiting(bool stretch = true, bool shear = true, bool bending = false);
	void setParameters_StrainLimiting(real min, real max, size_t maxIter);

	void postProcess(real dt, TriDynamicMesh<real>&) const;
	void postProcess(real dt, QuadDynamicMesh<real>&) const;

	const real strainLimiting(Vec3Array<real>& pos, const Constraints<real>&, const vector<Array<Spring<real>>*>& ss) const;

private : 
	bool		modifyStretch;
	bool		modifyShear;
	bool		modifyBending;

	real		minlimit, maxlimit;
	size_t		maxIter;
};

class NATUREDLL BridsonLimiting : public PostProcessingDynamics {
public : 
	BridsonLimiting(bool stretch = true, bool shear = true, bool bending = false);

	void setParameters_StrainLimiting(real min, real max, size_t maxIter);
	void setParameters_StrainRateLimiting(real limit, size_t maxIter);

	void postProcess(real dt, TriDynamicMesh<real>&) const;
	void postProcess(real dt, QuadDynamicMesh<real>&) const;

	const real strainLimiting(Vec3Array<real>& p, Vec3Array<real>& v, real dt, const Constraints<real>&, const vector<Array<Spring<real>>*>& ss) const;
	const real strainRateLimiting(Vec3Array<real>& p, Vec3Array<real>& v, real dt, const Constraints<real>&, const vector<Array<Spring<real>>*>& ss) const;

private : 
	bool		modifyStretch;
	bool		modifyShear;
	bool		modifyBending;

	real		minlimitStrain, maxlimitStrain;
	real		limitStrainRate;

	size_t		maxIter_Strain;
	size_t		maxIter_StrainRate;
};

};		// springspace

#endif