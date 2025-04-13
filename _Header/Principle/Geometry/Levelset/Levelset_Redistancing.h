#ifndef _GEOMETRYSPACE_LEVELSET_REDISTANCING_H_
#define _GEOMETRYSPACE_LEVELSET_REDISTANCING_H_

#include "Levelset2D.h"
#include "Levelset3D.h"
#include "FMMHeap.h"
#include "UpwindDerivative_Functor.h"
#include "TVD_RungeKutta.h"

namespace geometryspace {

using arithmetic::grid_index_type;
using std::vector;

template<typename T>
class Levelset_Redistancing {
public : 
	typedef typename Grid2<T,T>::index_type index_type;

	Levelset_Redistancing(T band) : narrowBandWidth(band) {}

	virtual void post_process() const { log_info << "Levelset_Redistancing _ Have no Post Process" << std::endl; }

	virtual void redistance(Levelset2D<T>&) const { basis::throwError("Not Implemented _ Levelset_Redistancing::redistance(Levelset2D<T>&)"); }	
	virtual void redistance(Levelset3D<T>&) const { basis::throwError("Not Implemented _ Levelset_Redistancing::redistance(Levelset3D<T>&)"); }	

	virtual void redistance(Grid2<T,T>&) const { basis::throwError("Not Implemented _ Levelset_Redistancing::redistance(Grid2<T>&)"); }	
	virtual void redistance(Grid3<T,T>&) const { basis::throwError("Not Implemented _ Levelset_Redistancing::redistance(Grid3<T>&)"); }	

	void setPhi_asExtreme(Grid2<T,T>&) const;
	void setPhi_asExtreme(Grid3<T,T>&) const;

	T	narrowBandWidth;			// Normalized Value	(0,1,2,...)		

	virtual void print_log(Log& log_) const { log_ << XmlElementName().c_str() << " _ Narrowband : " << narrowBandWidth << std::endl; }
	static const std::string XmlElementName() { return std::string("Levelset_Redistancing"); }
};

template<typename T>
class Levelset_RedistancingFMM : public Levelset_Redistancing<T> {
public : 
	Levelset_RedistancingFMM(T band) : Levelset_Redistancing(band), heap(NULL) {}
	virtual ~Levelset_RedistancingFMM() { if(heap) delete heap; }

	void redistance(Levelset2D<T>&) const;	
	void redistance(Levelset3D<T>&) const;	

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; Levelset_Redistancing::print_log(log_); }
	static const std::string XmlElementName() { return std::string("Levelset_RedistancingFMM"); }
	
private : 
	void redistance_initialize(Levelset2D<T>&) const;
	void redistance_initialize(Levelset3D<T>&) const;

	const T recomputeDistance(const Levelset2D<T>& lv, grid_index_type i, grid_index_type j) const;	
	const T recomputeDistance(Levelset3D<T>& lv, grid_index_type i, grid_index_type j, grid_index_type k) const;

	//////////////////////////////////////////////////////////////////////////
	// Variable for Fast Marching
	mutable FMMHeap<T>*			heap;
	mutable grid_index_type		minHeapSize;

	mutable Grid2<char,T>		sTag2D;	
	mutable Grid3<char,T>		sTag3D;		

	mutable vector<grid_index_type>	interfaceCellList;	
};

template<typename T, typename RK = RungeKutta_Euler>
class Levelset_Redistancing_PDE : public Levelset_Redistancing<T> {
public : 
	Levelset_Redistancing_PDE(T band, T _cfl, std::size_t iter = 0, bool peng = true) 
	: Levelset_Redistancing(band), cfl(_cfl), min_iter(iter), use_sign_peng(peng) { assert(cfl <= 0.5f); cfl = MIN(cfl, T(0.5)); }

	void set_min_iteration(std::size_t iter) { min_iter = iter; }

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << " _ CFL : " << cfl << " _ Peng Sign : " << use_sign_peng << std::endl; Levelset_Redistancing::print_log(log_); }
	static const std::string XmlElementName() { return std::string("Levelset_RedistancingPDE"); }

protected :
	static const T sign_sussman94(const T phi, const T h_square);

	static const T sign_peng99(const Grid2<T,T>& phi, const index_type i, const index_type j);
	static const T sign_peng99(const Grid3<T,T>& phi, const index_type i, const index_type j, const index_type k);

	template<typename SWEEP> void redistance_pde(Grid2<T,T>&, const Grid2<T,T>&, const SWEEP&) const;	
	template<typename SWEEP> void redistance_pde(Grid3<T,T>&, const Grid3<T,T>&, const SWEEP&) const;

	template<typename DEV> void redistance_pde_one_sweep(Grid2<T,T>&, const Grid2<T,T>&, const T tau) const;	
	template<typename DEV> void redistance_pde_one_sweep(Grid3<T,T>&, const Grid3<T,T>&, const T tau) const;

	T				cfl;
	std::size_t		min_iter;

	bool			use_sign_peng;

protected :
	template<typename DEV>
	class Redistancing_EE {		
	public :
		typedef DEV							dev_type;
		typedef typename DEV::value_type	value_type;

		Redistancing_EE(const Levelset_Redistancing_PDE& _redist) : redist(_redist), tau(0) {}

		void set_tau(value_type t) const { tau = t; }
		void operator()(typename DEV::grid_type& s, const typename DEV::grid_type& s0) const { redist.redistance_pde_one_sweep<DEV>(s, s0, tau); }

		const Levelset_Redistancing_PDE&	redist;		
		mutable value_type					tau;		
	};
};

template<typename T, typename RK = RungeKutta_Euler>
class Levelset_Redistancing_PDE_1stUpwind : public Levelset_Redistancing_PDE<T,RK> {
public : 
	Levelset_Redistancing_PDE_1stUpwind(T band, T _cfl, std::size_t min_iter = 0, bool _peng = true) : Levelset_Redistancing_PDE(band, _cfl, min_iter, _peng) {}

	void redistance(Levelset2D<T>& lv) const { redistance_pde(lv.phi, lv.phi, Redistancing_EE<Upwind_1st<Grid2<T,T>>>(*this)); }	
	void redistance(Levelset3D<T>& lv) const { redistance_pde(lv.phi, lv.phi, Redistancing_EE<Upwind_1st<Grid3<T,T>>>(*this)); }

	void redistance(Grid2<T,T>& phi) const { redistance_pde(phi, phi, Redistancing_EE<Upwind_1st<Grid2<T,T>>>(*this)); }
	void redistance(Grid3<T,T>& phi) const { redistance_pde(phi, phi, Redistancing_EE<Upwind_1st<Grid3<T,T>>>(*this)); }

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; Levelset_Redistancing_PDE::print_log(log_); }
	static const std::string XmlElementName() { return std::string("Levelset_Redistancing_1stUpwind"); }
};

template<typename T, typename RK = RungeKutta_Euler>
class Levelset_Redistancing_PDE_2ndUpwind : public Levelset_Redistancing_PDE<T,RK> {
public : 
	Levelset_Redistancing_PDE_2ndUpwind(T band, T _cfl, std::size_t min_iter = 0, bool _peng = true) : Levelset_Redistancing_PDE(band, _cfl, min_iter, _peng) {}

	void redistance(Levelset2D<T>& lv) const { redistance_pde(lv.phi, lv.phi, Redistancing_EE<Upwind_2nd<Grid2<T,T>>>(*this)); } 
	void redistance(Levelset3D<T>& lv) const { redistance_pde(lv.phi, lv.phi, Redistancing_EE<Upwind_2nd<Grid3<T,T>>>(*this)); } 

	void redistance(Grid2<T,T>& phi) const { redistance_pde(phi, phi, Redistancing_EE<Upwind_2nd<Grid2<T,T>>>(*this)); }
	void redistance(Grid3<T,T>& phi) const { redistance_pde(phi, phi, Redistancing_EE<Upwind_2nd<Grid3<T,T>>>(*this)); }

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; Levelset_Redistancing_PDE::print_log(log_); }

	static const std::string XmlElementName() { return std::string("Levelset_Redistancing_2ndUpwind"); }
};

template<typename T, typename RK = RungeKutta_Euler>
class Levelset_Redistancing_PDE_ENO : public Levelset_Redistancing_PDE<T,RK> {
public : 
	Levelset_Redistancing_PDE_ENO(T band, T _cfl, std::size_t min_iter = 0, bool _peng = true) : Levelset_Redistancing_PDE(band, _cfl, min_iter, _peng) {}

	void redistance(Levelset2D<T>& lv) const { redistance_pde(lv.phi, lv.phi, Redistancing_EE<Upwind_ENO<Grid2<T,T>>>(*this)); } 
	void redistance(Levelset3D<T>& lv) const { redistance_pde(lv.phi, lv.phi, Redistancing_EE<Upwind_ENO<Grid3<T,T>>>(*this)); } 

	void redistance(Grid2<T,T>& phi) const { redistance_pde(phi, phi, Redistancing_EE<Upwind_ENO<Grid2<T,T>>>(*this)); }
	void redistance(Grid3<T,T>& phi) const { redistance_pde(phi, phi, Redistancing_EE<Upwind_ENO<Grid3<T,T>>>(*this)); }

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; Levelset_Redistancing_PDE::print_log(log_); }
	static const std::string XmlElementName() { return std::string("Levelset_Redistancing_ENO"); }
};

template<typename T, typename RK = RungeKutta_Euler>
class Levelset_Redistancing_PDE_WENO : public Levelset_Redistancing_PDE<T,RK> {
public : 
	Levelset_Redistancing_PDE_WENO(T band, T _cfl, std::size_t min_iter = 0, bool _peng = true) : Levelset_Redistancing_PDE(band, _cfl, min_iter, _peng) {}

	void redistance(Levelset2D<T>& lv) const { redistance_pde(lv.phi, lv.phi, Redistancing_EE<Upwind_WENO<Grid2<T,T>>>(*this)); }
	void redistance(Levelset3D<T>& lv) const { redistance_pde(lv.phi, lv.phi, Redistancing_EE<Upwind_WENO<Grid3<T,T>>>(*this)); }

	void redistance(Grid2<T,T>& phi) const { redistance_pde(phi, phi, Redistancing_EE<Upwind_WENO<Grid2<T,T>>>(*this)); }
	void redistance(Grid3<T,T>& phi) const { redistance_pde(phi, phi, Redistancing_EE<Upwind_WENO<Grid3<T,T>>>(*this)); }

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; Levelset_Redistancing_PDE::print_log(log_); }
	static const std::string XmlElementName() { return std::string("Levelset_Redistancing_WENO"); }
};

};		// geometryspace

#include "Levelset_Redistancing_Utility.h"
#include "Levelset_RedistancingFMM2D.h"
#include "Levelset_RedistancingFMM3D.h"
#include "Levelset_Redistancing_PDE.h"

#endif