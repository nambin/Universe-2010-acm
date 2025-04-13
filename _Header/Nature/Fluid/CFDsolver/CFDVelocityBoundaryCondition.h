#ifndef _FLUIDSPACE_CFD_VELOCITY_BOUNDARY_CONDITION_H_
#define _FLUIDSPACE_CFD_VELOCITY_BOUNDARY_CONDITION_H_

#include "../../Nature.h"
#include "../../../Basis/ParallelProgramming/ParallelComputable.h"
#include "../../../Principle/Arithmetic/GridLibrary/GridBoundaryCondition.h"
#include "../Environment/Environments.h"
#include "../CFDsolver/CFDsolver2D.h"
#include "../CFDsolver/CFDsolver3D.h"

namespace fluidspace {

using arithmetic::GridBoundaryCondition;
using arithmetic::GridBoundaryCondition_Neumann;

class NATUREDLL CFDVelocityBoundaryCondition : public basis::ParallelComputable {
public : 
	CFDVelocityBoundaryCondition() : uBdry(NULL), vBdry(NULL), wBdry(NULL) {}

	virtual void post_process() const { log_info << "CFDVelocityBoundaryCondition _ Have no Post Process" << std::endl; }

	void setComputingOption(const ParallelComputable& rhs);

	void enforce(CFDsolver2D&);
	void enforce(CFDsolver3D&);

	GridBoundaryCondition<real,real>*	uBdry;
	GridBoundaryCondition<real,real>*	vBdry;
	GridBoundaryCondition<real,real>*	wBdry;

	void print_log(Log& log_) const { 
		if(uBdry) uBdry->print_log(log_);
		if(vBdry) vBdry->print_log(log_);
		if(wBdry) wBdry->print_log(log_);
	}
};

class NATUREDLL GridBoundaryCondition_U_NoSlip : public GridBoundaryCondition<real,real> {
public : 
	GridBoundaryCondition_U_NoSlip() {}
	void enforce(Grid2<real,real>&) const;
	void enforce(Grid3<real,real>&) const;

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; }

	static const std::string XmlElementName() { return std::string("GridBoundaryCondition_U_NoSlip"); }
};

class NATUREDLL GridBoundaryCondition_V_NoSlip : public GridBoundaryCondition<real,real> {
public : 
	GridBoundaryCondition_V_NoSlip() {}
	void enforce(Grid2<real,real>&) const;
	void enforce(Grid3<real,real>&) const;

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; }

	static const std::string XmlElementName() { return std::string("GridBoundaryCondition_V_NoSlip"); }
};

class NATUREDLL GridBoundaryCondition_W_NoSlip : public GridBoundaryCondition<real,real> {
public : 
	GridBoundaryCondition_W_NoSlip() {}
	void enforce(Grid2<real,real>&) const { throwError("Error in GridBoundaryCondition_W_NoSlip::enforce"); }
	void enforce(Grid3<real,real>&) const;

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; }

	static const std::string XmlElementName() { return std::string("GridBoundaryCondition_W_NoSlip"); }
};

class NATUREDLL GridBoundaryCondition_U_FreeSlip : public GridBoundaryCondition<real,real> {
public : 
	GridBoundaryCondition_U_FreeSlip() {}
	void enforce(Grid2<real,real>&) const;
	void enforce(Grid3<real,real>&) const;

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; }

	static const std::string XmlElementName() { return std::string("GridBoundaryCondition_U_FreeSlip"); }
};

class NATUREDLL GridBoundaryCondition_V_FreeSlip : public GridBoundaryCondition<real,real> {
public : 
	GridBoundaryCondition_V_FreeSlip() {}
	void enforce(Grid2<real,real>&) const;
	void enforce(Grid3<real,real>&) const;

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; }

	static const std::string XmlElementName() { return std::string("GridBoundaryCondition_V_FreeSlip"); }
};

class NATUREDLL GridBoundaryCondition_W_FreeSlip : public GridBoundaryCondition<real,real> {
public : 
	GridBoundaryCondition_W_FreeSlip() {}
	void enforce(Grid2<real,real>&) const { throwError("Error in GridBoundaryCondition_W_FreeSlip::enforce"); }
	void enforce(Grid3<real,real>&) const;

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; }

	static const std::string XmlElementName() { return std::string("GridBoundaryCondition_W_FreeSlip"); }
};

class NATUREDLL GridBoundaryCondition_V_CeilingSeparation : public GridBoundaryCondition<real,real> {
public : 
	GridBoundaryCondition_V_CeilingSeparation(real r) { ratio = r; }
	void enforce(Grid2<real,real>&) const;
	void enforce(Grid3<real,real>&) const;

	real ratio;

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << " _ Ratio : " << ratio << endl; }
	static const std::string XmlElementName() { return std::string("GridBoundaryCondition_V_CeilingSeparation"); }
};

};		// fluidspace

#endif