#ifndef _FLUIDSPACE_CFD_ADVECTION_H_
#define _FLUIDSPACE_CFD_ADVECTION_H_

#include "../../../Basis/Basis/Environments.h"
#include "../../../Principle/Arithmetic/GridLibrary/GridBoundaryCondition.h"
#include "../../../Principle/Geometry/MacGrid/MacGridFunctor.h"
#include "../../Nature.h"
#include "../Environment/Environments.h"
#include "../CFDsolver/CFDsolver2D.h"
#include "../CFDsolver/CFDsolver3D.h"
#include "../CFDsolver/CFDVelocityBoundaryCondition.h"
#include "../CFDsolver/CFDTracking.h"
#include "../DensityTracking/DensityTracking.h"
#include "../InterfaceTracking/InterfaceTracking_Levelset2D.h"
#include "../InterfaceTracking/InterfaceTracking_Levelset3D.h"
#include "../InterfaceTracking/InterfaceTracking_SRL2D.h"
#include "../InterfaceTracking/InterfaceTracking_SRL3D.h"
#include "../FluidLibrary/VelocityFields.h"

namespace fluidspace {

using arithmetic::GridBoundaryCondition;
using geometryspace::MacGrid2D;
using geometryspace::Levelset2D;
using geometryspace::Levelset3D;

class NATUREDLL CFDAdvection {
public : 
	typedef CFDsolver2D::index_type index_type;

public : 
	CFDAdvection() {}

	virtual void post_process() const { log_info << "CFDAdvection _ Have no Post Process" << std::endl; }

	virtual void advect(CFDsolver2D&, real dt, const CFDVelocityBoundaryCondition&) const;
	virtual void advect(CFDsolver3D&, real dt, const CFDVelocityBoundaryCondition&) const;

	virtual void advect(DensityTracking2D&, const CFDsolver2D&, real dt) const;	
	virtual void advect(DensityTracking3D&, const CFDsolver3D&, real dt) const;

	virtual void advect(InterfaceTracking_Levelset2D&, const CFDsolver2D&, real dt) const;	
	virtual void advect(InterfaceTracking_Levelset3D&, const CFDsolver3D&, real dt) const;

	virtual void advect(InterfaceTracking_Levelset2D&, const VelocityField2D<real>&, real dt) const;	
	virtual void advect(InterfaceTracking_Levelset3D&, const VelocityField3D<real>&, real dt) const;

	virtual void advect(InterfaceTracking_SRL2D&, const CFDsolver2D&, real dt) const { basis::throwError("Not Implemented _ advect SRL"); }
	virtual void advect(InterfaceTracking_SRL2D&, const VelocityField2D<real>&, real dt) const { basis::throwError("Not Implemented _ advect SRL"); }

	virtual void advect(InterfaceTracking_SRL3D&, const CFDsolver3D&, real dt) const { basis::throwError("Not Implemented _ advect SRL"); }
	virtual void advect(InterfaceTracking_SRL3D&, const VelocityField3D<real>&, real dt) const { basis::throwError("Not Implemented _ advect SRL"); }

	virtual void advect(Grid2r& s, const Grid2r& s0, const CFDsolver2D&, real dt, MacGrid2D<real>::FIELD_TYPE type, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const = 0;
	virtual void advect(Grid3r& s, const Grid3r& s0, const CFDsolver3D&, real dt, MacGrid3D<real>::FIELD_TYPE type, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const = 0;

	virtual void advect(Grid2r& s, const Grid2r& s0, const VelocityField2D<real>&, real dt, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const = 0;
	virtual void advect(Grid3r& s, const Grid3r& s0, const VelocityField3D<real>&, real dt, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const = 0;

	virtual void print_log(Log& log_) const {}

	static const std::string XmlElementName() { return std::string("CFDAdvection"); }
};

template<typename TRACER = CFDTracking_EE>
class CFDAdvection_SemiLagrangian : public CFDAdvection {
public : 
	CFDAdvection_SemiLagrangian() {}

	void advect(Grid2r& s, const Grid2r& s0, const CFDsolver2D& cfd, real dt, MacGrid2D<real>::FIELD_TYPE type = MacGrid2D<real>::TYPE_S, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const { advect_sl(s,s0,cfd,dt,type,b); }
	void advect(Grid3r& s, const Grid3r& s0, const CFDsolver3D& cfd, real dt, MacGrid3D<real>::FIELD_TYPE type = MacGrid3D<real>::TYPE_S, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const { advect_sl(s,s0,cfd,dt,type,b); }

	void advect(Grid2r& s, const Grid2r& s0, const VelocityField2D<real>& vel, real dt, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const { advect_sl(s,s0,vel,dt,b); }
	void advect(Grid3r& s, const Grid3r& s0, const VelocityField3D<real>& vel, real dt, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const { advect_sl(s,s0,vel,dt,b); }

	static void advect_sl(Grid2r& s, const Grid2r& s0, const CFDsolver2D&, real dt, MacGrid2D<real>::FIELD_TYPE type, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>());	
	static void advect_sl(Grid3r& s, const Grid3r& s0, const CFDsolver3D&, real dt, MacGrid3D<real>::FIELD_TYPE type, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>());

	static void advect_sl(Grid2r& s, const Grid2r& s0, const VelocityField2D<real>&, real dt, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>());	
	static void advect_sl(Grid3r& s, const Grid3r& s0, const VelocityField3D<real>&, real dt, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>());

	static void limiting_sl(Grid2r& s, const Grid2r& s0, const CFDsolver2D&, real dt, MacGrid2D<real>::FIELD_TYPE type);
	static void limiting_sl(Grid3r& s, const Grid3r& s0, const CFDsolver3D&, real dt, MacGrid3D<real>::FIELD_TYPE type);

	static void limiting_sl(Grid2r& s, const Grid2r& s0, const VelocityField2D<real>&, real dt);		
	static void limiting_sl(Grid3r& s, const Grid3r& s0, const VelocityField3D<real>&, real dt);	

	// Revert to Semi-Lagrangian (if extrema is occurred)	
	// s is result of Semi-Lagrangian
	// sTmp is result of other Solver (ex. BFECC, MacCormack)
	static void revert_sl(Grid2r& s, const Grid2r& s0, const Grid2r& sTmp, const CFDsolver2D&, real dt, MacGrid2D<real>::FIELD_TYPE type);
	static void revert_sl(Grid2r& s, const Grid2r& s0, const Grid2r& sTmp, const VelocityField2D<real>&, real dt);	

	static void revert_sl(Grid3r& s, const Grid3r& s0, const Grid3r& sTmp, const CFDsolver3D&, real dt, MacGrid3D<real>::FIELD_TYPE type);
	static void revert_sl(Grid3r& s, const Grid3r& s0, const Grid3r& sTmp, const VelocityField3D<real>&, real dt);	

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; CFDAdvection::print_log(log_); }
	static const std::string XmlElementName() { return std::string("CFDAdvection_SemiLagrangian"); }
};

template<typename TRACER = CFDTracking_RK2>
class CFDAdvection_SemiLagrangian_SRL : public CFDAdvection_SemiLagrangian<TRACER> {
public : 
	CFDAdvection_SemiLagrangian_SRL(bool clip = true) { clipping = clip; time_total = time_only_advection = 0.0; num_itrp = num_clipping = 0; }

	void post_process() const;

	void advect(InterfaceTracking_SRL2D& s, const CFDsolver2D& cfd, real dt) const { advect_srl(s, cfd, dt); }
	void advect(InterfaceTracking_SRL2D& s, const VelocityField2D<real>& vel, real dt) const { advect_srl(s, vel, dt); }

	virtual void advect(InterfaceTracking_SRL3D& s, const CFDsolver3D& cfd, real dt) const { advect_srl(s, cfd, dt); }
	virtual void advect(InterfaceTracking_SRL3D& s, const VelocityField3D<real>& vel, real dt) const { advect_srl(s, vel, dt); }

	template<typename CFD> void advect_srl(InterfaceTracking_SRL2D&, const CFD&, real dt) const;
	template<typename CFD> void advect_srl(InterfaceTracking_SRL3D&, const CFD&, real dt) const;

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; CFDAdvection::print_log(log_); }
	static const std::string XmlElementName() { return std::string("CFDAdvection_SemiLagrangian_SRL"); }

	bool clipping;

	mutable double			time_total, time_only_advection;
	mutable unsigned int	num_itrp, num_clipping;
};

class NATUREDLL CFDAdvection_BFECC : public CFDAdvection {
public : 
	CFDAdvection_BFECC(bool _limit = false) : limiting(_limit) {}

	void advect(Grid2r& s, const Grid2r& s0, const CFDsolver2D&, real dt, MacGrid2D<real>::FIELD_TYPE type, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;	
	void advect(Grid3r& s, const Grid3r& s0, const CFDsolver3D&, real dt, MacGrid3D<real>::FIELD_TYPE type, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;

	void advect(Grid2r& s, const Grid2r& s0, const VelocityField2D<real>&, real dt, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;	
	void advect(Grid3r& s, const Grid3r& s0, const VelocityField3D<real>&, real dt, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; CFDAdvection::print_log(log_); }
	static const std::string XmlElementName() { return std::string("CFDAdvection_BFECC"); }

	bool limiting;
};

class NATUREDLL CFDAdvection_MacCormack : public CFDAdvection {
public : 
	CFDAdvection_MacCormack(bool _limit = true, bool _revertsl = true) : limiting(_limit), revertsl(_revertsl) {}

	void advect(Grid2r& s, const Grid2r& s0, const CFDsolver2D&, real dt, MacGrid2D<real>::FIELD_TYPE type, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;	
	void advect(Grid3r& s, const Grid3r& s0, const CFDsolver3D&, real dt, MacGrid3D<real>::FIELD_TYPE type, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;

	void advect(Grid2r& s, const Grid2r& s0, const VelocityField2D<real>&, real dt, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;	
	void advect(Grid3r& s, const Grid3r& s0, const VelocityField3D<real>&, real dt, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;

	bool limiting, revertsl;
};

class CFDAdvection_Upwind_Conservative : public CFDAdvection {
public : 
	CFDAdvection_Upwind_Conservative(real _cfl) { cfl = _cfl; }

	void advect(Grid2r& s, const Grid2r& s0, const CFDsolver2D&, real dt, MacGrid2D<real>::FIELD_TYPE type, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;	
	void advect(Grid3r& s, const Grid3r& s0, const CFDsolver3D&, real dt, MacGrid3D<real>::FIELD_TYPE type, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;

	void advect(Grid2r& s, const Grid2r& s0, const VelocityField2D<real>&, real dt, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;	
	void advect(Grid3r& s, const Grid3r& s0, const VelocityField3D<real>&, real dt, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;

	template<class CFD> void advect_ee(Grid2r& s, const Grid2r& s0, const CFD&, real dt, const GridBdry& b) const;	
	template<class CFD> void advect_ee(Grid3r& s, const Grid3r& s0, const CFD&, real dt, const GridBdry& b) const;	

	void print_log(Log& log_) const { log_ << XmlElementName().c_str(); log_ << " _ CFL : " << cfl << endl; CFDAdvection::print_log(log_); }
	static const std::string XmlElementName() { return std::string("CFDAdvection_Upwind_Conservative"); }

	real cfl;
};

template<typename RK>
class CFDAdvection_Upwind : public CFDAdvection {
public : 
	CFDAdvection_Upwind(real _cfl) { cfl = _cfl; }

	template<typename DEV> static void advect_ee(Grid2r& s, const Grid2r& s0, const CFDsolver2D&, real dt, MacGrid2D<real>::FIELD_TYPE type, const GridBdry& b, real cfl);	
	template<typename DEV> static void advect_ee(Grid3r& s, const Grid3r& s0, const CFDsolver3D&, real dt, MacGrid3D<real>::FIELD_TYPE type, const GridBdry& b, real cfl);

	template<typename DEV> static void advect_ee(Grid2r& s, const Grid2r& s0, const VelocityField2D<real>&, real dt, const GridBdry& b, real cfl);	
	template<typename DEV> static void advect_ee(Grid3r& s, const Grid3r& s0, const VelocityField3D<real>&, real dt, const GridBdry& b, real cfl);

	real cfl;

protected :
	template<typename VEL, typename DEV>
	class Advect_EE_VelField {
	public :
		Advect_EE_VelField(const VEL& _vel, real _dt, const GridBdry& _b, real _cfl) : vel(_vel), dt(_dt), b(_b), cfl(_cfl) {}

		void operator()(typename DEV::grid_type& s, const typename DEV::grid_type& s0) const { CFDAdvection_Upwind::advect_ee<DEV>(s, s0, vel, dt, b, cfl);		}

		const VEL&					vel;
		real						dt, cfl;
		const GridBdry&				b;
	};

	template<typename CFD, typename DEV>
	class Advect_EE_CFD {
	public :
		Advect_EE_CFD(const CFD& _cfd, real _dt, typename CFD::FIELD_TYPE _type, const GridBdry& _b, real _cfl) : cfd(_cfd), type(_type), dt(_dt), b(_b), cfl(_cfl) {}

		void operator()(typename DEV::grid_type& s, const typename DEV::grid_type& s0) const { CFDAdvection_Upwind::advect_ee<DEV>(s, s0, cfd, dt, type, b, cfl); }

		const CFD&					cfd;
		typename CFD::FIELD_TYPE	type;
		real						dt, cfl;
		const GridBdry&				b;
	};
};

template<typename RK>
class CFDAdvection_1stUpwind : public CFDAdvection_Upwind<RK> {
public : 
	CFDAdvection_1stUpwind(real cfl = 1.0e+6) : CFDAdvection_Upwind(cfl) {}

	void advect(Grid2r& s, const Grid2r& s0, const CFDsolver2D&, real dt, MacGrid2D<real>::FIELD_TYPE type, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;	
	void advect(Grid3r& s, const Grid3r& s0, const CFDsolver3D&, real dt, MacGrid3D<real>::FIELD_TYPE type, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;

	void advect(Grid2r& s, const Grid2r& s0, const VelocityField2D<real>&, real dt, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;	
	void advect(Grid3r& s, const Grid3r& s0, const VelocityField3D<real>&, real dt, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;

	void print_log(Log& log_) const { log_ << XmlElementName().c_str(); log_ << " _ CFL : " << cfl << endl; CFDAdvection::print_log(log_); }
	static const std::string XmlElementName() { return std::string("CFDAdvection_1stUpwind"); }
};

template<typename RK>
class CFDAdvection_2ndUpwind : public CFDAdvection_Upwind<RK> {
public : 
	CFDAdvection_2ndUpwind(real cfl = 1.0e+6) : CFDAdvection_Upwind(cfl) {}

	void advect(Grid2r& s, const Grid2r& s0, const CFDsolver2D&, real dt, MacGrid2D<real>::FIELD_TYPE type, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;	
	void advect(Grid3r& s, const Grid3r& s0, const CFDsolver3D&, real dt, MacGrid3D<real>::FIELD_TYPE type, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;

	void advect(Grid2r& s, const Grid2r& s0, const VelocityField2D<real>&, real dt, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;	
	void advect(Grid3r& s, const Grid3r& s0, const VelocityField3D<real>&, real dt, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;

	void print_log(Log& log_) const { log_ << XmlElementName().c_str(); log_ << " _ CFL : " << cfl << endl; CFDAdvection::print_log(log_); }
	static const std::string XmlElementName() { return std::string("CFDAdvection_2ndUpwind"); }
};

template<typename RK>
class CFDAdvection_ENO : public CFDAdvection_Upwind<RK> {
public : 
	CFDAdvection_ENO(real cfl = 1.0e+6) : CFDAdvection_Upwind(cfl) {}

	void advect(Grid2r& s, const Grid2r& s0, const CFDsolver2D&, real dt, MacGrid2D<real>::FIELD_TYPE type, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;	
	void advect(Grid3r& s, const Grid3r& s0, const CFDsolver3D&, real dt, MacGrid3D<real>::FIELD_TYPE type, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;

	void advect(Grid2r& s, const Grid2r& s0, const VelocityField2D<real>&, real dt, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;	
	void advect(Grid3r& s, const Grid3r& s0, const VelocityField3D<real>&, real dt, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;

	void print_log(Log& log_) const { log_ << XmlElementName().c_str(); log_ << " _ CFL : " << cfl << endl; CFDAdvection::print_log(log_); }
	static const std::string XmlElementName() { return std::string("CFDAdvection_ENO"); }
};

template<typename RK>
class CFDAdvection_WENO : public CFDAdvection_Upwind<RK> {
public : 
	CFDAdvection_WENO(real cfl = 1.0e+6) : CFDAdvection_Upwind(cfl) {}

	void advect(Grid2r& s, const Grid2r& s0, const CFDsolver2D&, real dt, MacGrid2D<real>::FIELD_TYPE type, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;	
	void advect(Grid3r& s, const Grid3r& s0, const CFDsolver3D&, real dt, MacGrid3D<real>::FIELD_TYPE type, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;

	void advect(Grid2r& s, const Grid2r& s0, const VelocityField2D<real>&, real dt, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;	
	void advect(Grid3r& s, const Grid3r& s0, const VelocityField3D<real>&, real dt, const GridBdry& b = GridBoundaryCondition_Neumann<real,real>()) const;

	void print_log(Log& log_) const { log_ << XmlElementName().c_str(); log_ << " _ CFL : " << cfl << endl; CFDAdvection::print_log(log_); }
	static const std::string XmlElementName() { return std::string("CFDAdvection_WENO"); }
};

//class NATUREDLL CFDAdvection_USCIP : public CFDAdvection {
//public : 
//	CFDAdvection_USCIP(bool _limit = false) : limiting(_limit) {}
//
//	void advect(Grid2r& s, Grid2r& sx, Grid2r& sy, const Grid2r& s0, const Grid2r& sx0, const Grid2r& sy0, const CFDsolver2D&, real dt, MacGrid2D<real>::FIELD_TYPE type);	
//	void advect(Grid3r& s, Grid3r& sx, Grid3r& sy, const Grid3r& s0, const Grid3r& sx0, const Grid3r& sy0, const CFDsolver3D&, real dt, MacGrid3D<real>::FIELD_TYPE type);
//
//	void advect(Grid2r& s, Grid2r& sx, Grid2r& sy, const Grid2r& s0, const Grid2r& sx0, const Grid2r& sy0, const VelocityField2D<real>&, real dt);	
//	void advect(Grid3r& s, Grid3r& sx, Grid3r& sy, const Grid3r& s0, const Grid3r& sx0, const Grid3r& sy0, const VelocityField3D<real>&, real dt);
//
//	bool limiting;
//};

};		// fluidspace

#include "CFDAdvection_SemiLagrangian.h"
#include "CFDAdvection_SRL.h"
#include "CFDAdvection_Upwind.h"
#include "CFDAdvection_Upwind_Conservative.h"

#endif