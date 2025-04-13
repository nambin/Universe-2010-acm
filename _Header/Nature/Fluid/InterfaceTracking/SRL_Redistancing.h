#pragma once

#include "../../../Principle/Geometry/Levelset/TVD_RungeKutta.h"
#include "../../../Principle/Geometry/Levelset/Levelset_Redistancing.h"
#include "../../../Principle/Geometry/Levelset/UpwindDerivative_Functor.h"
#include "../../../Principle/Geometry/SRL/SRL_UpwindDerivative_Functor.h"
#include "../../Nature.h"
#include "InterfaceTracking_SRL2D.h"
#include "InterfaceTracking_SRL3D.h"

namespace fluidspace {

class NATUREDLL SRL_Redistancing {
public : 
	typedef geometryspace::SRL<real>::index_type index_type;

	SRL_Redistancing() { cfl_frequency = 1; min_updates = 1; max_updates = 1; cfl_acc = 0; num_updates = 0; }

	virtual void post_process() const { log_info << "SRL_Redistancing _ Have no Post Process" << std::endl; }

	void set_redistance_cfl_frequency(real freq) { cfl_frequency = freq; }
	void set_redistance_min_max_updates(unsigned int _min, unsigned int _max) { min_updates = _min; max_updates = _max; }

	const bool is_times_up(real dt, real C_CFL) const;

	virtual void redistance(InterfaceTracking_SRL2D&) const { basis::throwError("Not Implemented _ SRL_Redistancing::redistance(SRL2D<T>&)"); }	
	virtual void redistance(InterfaceTracking_SRL3D&) const { basis::throwError("Not Implemented _ SRL_Redistancing::redistance(SRL3D<T>&)"); }	

	void importBinary(std::ifstream& ifs);
	void exportBinary(std::ofstream& ofs) const;

	virtual void print_log(Log& log_) const { 
		log_ << XmlElementName().c_str() << " _ CFL Frequency : " << cfl_frequency << " _ "; 
		log_ << "Min, Max Updates : " << min_updates << ", " << max_updates << endl;  
	}
	static const std::string XmlElementName() { return std::string("SRL_Redistancing"); }

	unsigned int			min_updates, max_updates;
	real					cfl_frequency;

	mutable real			cfl_acc;
	mutable unsigned int	num_updates;
};

class NATUREDLL SRL_Redistancing_No : public SRL_Redistancing {
public : 
	typedef geometryspace::SRL<real>::index_type index_type;

	SRL_Redistancing_No() {}

	void redistance(InterfaceTracking_SRL2D&) const {}
	void redistance(InterfaceTracking_SRL3D&) const {}

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; SRL_Redistancing::print_log(log_); }
	static const std::string XmlElementName() { return std::string("SRL_Redistancing_No"); }
};

class NATUREDLL SRL_Redistancing_Direct : public SRL_Redistancing, public geometryspace::Levelset_Redistancing_PDE<real> {
public : 
	typedef geometryspace::SRL<real>::index_type index_type;

	SRL_Redistancing_Direct(real _band) : band(_band), geometryspace::Levelset_Redistancing_PDE<real>(_band, 0.5, 0, true) {}

	void redistance(InterfaceTracking_SRL2D&) const;
	void redistance(InterfaceTracking_SRL3D&) const;

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; SRL_Redistancing::print_log(log_); }
	static const std::string XmlElementName() { return std::string("SRL_Redistancing_Direct"); }

protected : 
	const real get_distance(const Vector2<real>&, const SRL2D<real>&, const vector<vector<Segment2D<real>>>&) const;

	real band;			// Normalized value 0,1,2,...
};

class NATUREDLL SRL_Redistancing_FMM : public SRL_Redistancing {
public : 
	SRL_Redistancing_FMM(real band) {}
};

class NATUREDLL SRL_Redistancing_PDE : public SRL_Redistancing, public geometryspace::Levelset_Redistancing_PDE<real> {
public : 
	typedef SRL_Redistancing::index_type index_type;

	void post_process() const;

	SRL_Redistancing_PDE(real band, real _cfl, bool only_narrow = true, std::size_t iter = 0) 
	: geometryspace::Levelset_Redistancing_PDE<real>(band, _cfl, iter, false), narrowband_calculation_only(only_narrow) { time_total = time_only_fine_grid = 0; }

	template<typename DEV_SRL, typename DEV_LV> void redistance_pde(InterfaceTracking_SRL2D&) const;
	template<typename DEV_SRL, typename DEV_LV> void redistance_pde(InterfaceTracking_SRL3D&) const;

	void print_log(Log& log_) const { SRL_Redistancing::print_log(log_); geometryspace::Levelset_Redistancing_PDE<real>::print_log(log_); }
	static const std::string XmlElementName() { return std::string("SRL_Redistancing_PDE"); }

private : 
	template<typename DEV_SRL, typename DEV_LV> void redistance_pde_wholedomain(InterfaceTracking_SRL2D&) const;
	template<typename DEV_SRL, typename DEV_LV> void redistance_pde_narrowband(InterfaceTracking_SRL2D&) const;

	template<typename DEV_SRL, typename DEV_LV> void redistance_pde_wholedomain(InterfaceTracking_SRL3D&) const { basis::throwError("Not Implemented"); }
	template<typename DEV_SRL, typename DEV_LV> void redistance_pde_narrowband(InterfaceTracking_SRL3D&) const;

protected : 
	bool narrowband_calculation_only;

	mutable double	time_total, time_only_fine_grid;
};

class NATUREDLL SRL_Redistancing_PDE_1stUpwind : public SRL_Redistancing_PDE {
public : 
	SRL_Redistancing_PDE_1stUpwind(real band, real _cfl, bool only_narrow = true, std::size_t min_iter = 0) : SRL_Redistancing_PDE(band, _cfl, only_narrow, min_iter) {}

	void redistance(InterfaceTracking_SRL2D& srl) const { redistance_pde<Upwind_SRL_1st<SRL2D<real>>, Upwind_1st<Grid2<real,real>>>(srl); }
	void redistance(InterfaceTracking_SRL3D& srl) const { redistance_pde<Upwind_SRL_1st<SRL3D<real>>, Upwind_1st<Grid3<real,real>>>(srl); }

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; SRL_Redistancing_PDE::print_log(log_); }
	static const std::string XmlElementName() { return std::string("SRL_Redistancing_PDE_1stUpwind"); }
};

class NATUREDLL SRL_Redistancing_PDE_2ndUpwind : public SRL_Redistancing_PDE {
public : 
	SRL_Redistancing_PDE_2ndUpwind(real band, real _cfl, bool only_narrow = true, std::size_t min_iter = 0) : SRL_Redistancing_PDE(band, _cfl, only_narrow, min_iter) {}

	void redistance(InterfaceTracking_SRL2D& srl) const { redistance_pde<Upwind_SRL_2nd<SRL2D<real>>, Upwind_1st<Grid2<real,real>>>(srl); }
	void redistance(InterfaceTracking_SRL3D& srl) const { redistance_pde<Upwind_SRL_2nd<SRL3D<real>>, Upwind_1st<Grid3<real,real>>>(srl); }

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; SRL_Redistancing_PDE::print_log(log_); }
	static const std::string XmlElementName() { return std::string("SRL_Redistancing_PDE_2ndUpwind"); }
};

class NATUREDLL SRL_Redistancing_PDE_WENO3rd : public SRL_Redistancing_PDE {
public : 
	SRL_Redistancing_PDE_WENO3rd(real band, real _cfl, bool only_narrow = true, std::size_t min_iter = 0) : SRL_Redistancing_PDE(band, _cfl, only_narrow, min_iter) {}

	void redistance(InterfaceTracking_SRL2D& srl) const { redistance_pde<Upwind_SRL_WENO3rd<SRL2D<real>>, Upwind_1st<Grid2<real,real>>>(srl); }
	void redistance(InterfaceTracking_SRL3D& srl) const { redistance_pde<Upwind_SRL_WENO3rd<SRL3D<real>>, Upwind_1st<Grid3<real,real>>>(srl); }

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; SRL_Redistancing_PDE::print_log(log_); }
	static const std::string XmlElementName() { return std::string("SRL_Redistancing_PDE_WENO3rd"); }
};

class NATUREDLL SRL_Redistancing_PDE_WENO : public SRL_Redistancing_PDE {
public : 
	SRL_Redistancing_PDE_WENO(real band, real _cfl, bool only_narrow = true, std::size_t min_iter = 0) : SRL_Redistancing_PDE(band, _cfl, only_narrow, min_iter) {}

	void redistance(InterfaceTracking_SRL2D& srl) const { redistance_pde<Upwind_SRL_WENO<SRL2D<real>>, Upwind_1st<Grid2<real,real>>>(srl); }
	void redistance(InterfaceTracking_SRL3D& srl) const { redistance_pde<Upwind_SRL_WENO<SRL3D<real>>, Upwind_1st<Grid3<real,real>>>(srl); }

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; SRL_Redistancing_PDE::print_log(log_); }
	static const std::string XmlElementName() { return std::string("SRL_Redistancing_PDE_WENO"); }
};

};			// fluidspace

#include "SRL_Redistancing_PDE.h"

