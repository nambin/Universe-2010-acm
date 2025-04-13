#ifndef _FLUIDSPACE_CFD_PROJECTION_H_
#define _FLUIDSPACE_CFD_PROJECTION_H_

#include "../Environment/Environments.h"
#include "../CFDsolver/CFDsolver2D.h"
#include "../InterfaceTracking/InterfaceTracking_Levelset2D.h"
#include "../../Nature.h"
#include "../../../Principle/Arithmetic/GridLibrary/Grid2_PoissonSolver_CG.h"
#include "../../../Principle/Arithmetic/GridLibrary/Grid3_PoissonSolver_CG.h"
#include "../../../Principle/Arithmetic/GridLibrary/Grid2_PoissonSolver_GMG.h"
#include "../../../Principle/Arithmetic/GridLibrary/Grid3_PoissonSolver_GMG.h"
#include "../../../Basis/ParallelProgramming/ParallelComputable.h"

namespace fluidspace {

using arithmetic::Grid2_PoissonSolver;
using arithmetic::Grid3_PoissonSolver;

class NATUREDLL CFDProjection : public basis::ParallelComputable {
public : 
	typedef CFDsolver2D::index_type index_type;

	CFDProjection() : sol2D(NULL), sol3D(NULL) { thres = 1.0e-3; }

	virtual void post_process() const { log_info << "CFDProjection _ Have no Post Process" << std::endl; }

	void set_threshold(const double _thres) { thres = _thres; }
	void set_divergence_threshold(const double _div_thres) { if(sol2D) thres = _div_thres * sol2D->x.dx; else if(sol3D) thres = _div_thres * sol3D->x.dx; else throwError("Error in CFDProjection"); }

	void setComputingOption(const ParallelComputable& rhs);
	
	void init_PoissonSolver_CG(const CFDsolver2D&);
	void init_PoissonSolver_CG(const CFDsolver3D&);

	void setPreconditioner_ICCG();
	void setPreconditioner_BlockRedBlackICCG(index_type nbX, index_type nbY);
	void setPreconditioner_BlockRedBlackICCG(index_type nbX, index_type nbY, index_type nbZ);

	void init_PoissonSolver_GMG(const CFDsolver2D&, unsigned int num_level = 12);
	void init_PoissonSolver_GMG(const CFDsolver3D&, unsigned int num_level = 12);
	
	void set_iteration_fine_grid_GMG(size_t iter);
	void set_iteration_number_GMG(size_t pre, size_t post);
	void set_parallel_threshold_GMG(index_type thres);

	virtual void set_initial_guess(Grid2_PoissonSolver<double>&) = 0;
	virtual void set_initial_guess(Grid3_PoissonSolver<double>&) = 0;
	
	virtual size_t project(CFDsolver2D&, real dt, size_t iter = 300) = 0;
	virtual size_t project(CFDsolver3D&, real dt, size_t iter = 300) = 0;

	virtual void importBinary(std::ifstream& ifs);
	virtual void exportBinary(std::ofstream& ofs) const;

	virtual void print_log(Log& log_) const { log_ << "Projection Threshold _ " << thres << std::endl; }
	static const std::string XmlElementName() { return std::string("CFDProjection"); }

protected : 
	void buildDivergence(Grid2d&, const CFDsolver2D&);
	void buildDivergence(Grid3d&, const CFDsolver3D&);

	void addPressureGradient(CFDsolver2D&, const Grid2d& pressure);
	void addPressureGradient(CFDsolver3D&, const Grid3d& pressure);

	void buildPoissonMatrix(Grid2x3d&, Grid2<bool,double>&, const Grid2d& u_part, const Grid2d& v_part);
	void buildPoissonMatrix(Grid3x4d&, Grid3<bool,double>&, const Grid3d& u_part, const Grid3d& v_part, const Grid3d& w_part);

	void buildPoissonMatrix(Grid2x3d&, Grid2<bool,double>&, const Grid2d& u_part, const Grid2d& v_part, const Grid2d& s_part);
	void buildPoissonMatrix(Grid3x4d&, Grid3<bool,double>&, const Grid3d& u_part, const Grid3d& v_part, const Grid3d& w_part, const Grid3d& s_part);

	void buildPoissonMatrix(Grid2_PoissonSolver_GMG<double>&, const Grid2d& u_part, const Grid2d& v_part);
	void buildPoissonMatrix(Grid3_PoissonSolver_GMG<double>&, const Grid3d& u_part, const Grid3d& v_part, const Grid3d& w_part);

	void buildPoissonMatrix(Grid2_PoissonSolver_GMG<double>&, const Grid2d& u_part, const Grid2d& v_part, const Grid2d& s_part);
	void buildPoissonMatrix(Grid3_PoissonSolver_GMG<double>&, const Grid3d& u_part, const Grid3d& v_part, const Grid3d& w_part, const Grid3d& s_part);

	void coarsening_PoissonMatrix(Grid2d& u_coarse, Grid2d& v_coarse, const Grid2d& u_fine, const Grid2d& v_fine);
	void coarsening_PoissonMatrix(Grid3d& u_coarse, Grid3d& v_coarse, Grid3d& w_coarse, const Grid3d& u_fine, const Grid3d& v_fine, const Grid3d& w_fine);

	void coarsening_PoissonMatrix(Grid2d& u_coarse, Grid2d& v_coarse, Grid2d& s_coarse, const Grid2d& u_fine, const Grid2d& v_fine, const Grid2d& s_fine);
	void coarsening_PoissonMatrix(Grid3d& u_coarse, Grid3d& v_coarse, Grid3d& w_coarse, Grid3d& s_coarse, const Grid3d& u_fine, const Grid3d& v_fine, const Grid3d& w_fine, const Grid3d& s_fine);

	Grid2_PoissonSolver<double>* sol2D;
	Grid3_PoissonSolver<double>* sol3D;

	double thres;
};

class NATUREDLL CFDProjection_SinglePhase : public CFDProjection {
public : 
	CFDProjection_SinglePhase() {}

	void set_initial_guess(Grid2_PoissonSolver<double>&) {}
	void set_initial_guess(Grid3_PoissonSolver<double>&) {}

	size_t project(CFDsolver2D& cfd, real dt, size_t iter) { return project_SinglePhase(cfd, *sol2D, iter, thres); }
	size_t project(CFDsolver3D& cfd, real dt, size_t iter) { return project_SinglePhase(cfd, *sol3D, iter, thres); }

	template<typename CFD, typename SOLVER> size_t project_SinglePhase(CFD&, SOLVER& sol, size_t iter, double thres);

	void buildPoissonMatrix(Grid2d& u_part, Grid2d& v_part, const CFDsolver2D&);
	void buildPoissonMatrix(Grid3d& u_part, Grid3d& v_part, Grid3d& w_part, const CFDsolver3D&);	

	void buildPoissonMatrix(Grid2_PoissonSolver<double>&, const CFDsolver2D&);
	void buildPoissonMatrix(Grid3_PoissonSolver<double>&, const CFDsolver3D&);	

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; CFDProjection::print_log(log_); }
	static const std::string XmlElementName() { return std::string("CFDProjection_SinglePhase"); }
};

class NATUREDLL CFDProjection_Liquid : public CFDProjection {
public : 
	CFDProjection_Liquid() : tracking(NULL), dt(0.03f), liquid_biased_filter(false) {}

	void set_InterfaceTracking(const InterfaceTracking* track) { tracking = track; }

	void set_liquid_biased_filter(const bool _flag) { liquid_biased_filter = _flag; }

	template<typename CFD, typename SOLVER> size_t project_Liquid(CFD&, SOLVER& sol, real dt, size_t iter, double thres);

	virtual void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; CFDProjection::print_log(log_); }
	static const std::string XmlElementName() { return std::string("CFDProjection_Liquid"); }

protected : 
	virtual void buildDivergence(Grid2d&, const Grid2r& phi, const CFDsolver2D&) = 0;	
	virtual void buildDivergence(Grid3d&, const Grid3r& phi, const CFDsolver3D&) = 0;

	virtual void buildPoissonMatrix(Grid2_PoissonSolver<double>&, const Grid2r& phi, const CFDsolver2D&) = 0;
	virtual void buildPoissonMatrix(Grid3_PoissonSolver<double>&, const Grid3r& phi, const CFDsolver3D&) = 0;	

	virtual void buildPoissonMatrix(Grid2d& u_part, Grid2d& v_part, Grid2d& s_part, const Grid2r& phi, const CFDsolver2D&) = 0;
	virtual void buildPoissonMatrix(Grid3d& u_part, Grid3d& v_part, Grid3d& w_part, Grid3d& s_part, const Grid3r& phi, const CFDsolver3D&) = 0;

	virtual void addPressureGradient(CFDsolver2D& cfd, const Grid2d& pressure, const Grid2r& phi) = 0;
	virtual void addPressureGradient(CFDsolver3D& cfd, const Grid3d& pressure, const Grid3r& phi) = 0;

	const InterfaceTracking* tracking;
	
	bool			liquid_biased_filter;
	mutable real	dt;
};

class NATUREDLL CFDProjection_FreeSurface : public CFDProjection_Liquid {
public : 
	CFDProjection_FreeSurface() {}

	void set_initial_guess(Grid2_PoissonSolver<double>&);
	void set_initial_guess(Grid3_PoissonSolver<double>&);

	size_t project(CFDsolver2D& cfd, real dt, size_t iter) { return project_Liquid(cfd, *sol2D, dt, iter, thres); }
	size_t project(CFDsolver3D& cfd, real dt, size_t iter) { return project_Liquid(cfd, *sol3D, dt, iter, thres); }

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << endl; CFDProjection_Liquid::print_log(log_); }
	static const std::string XmlElementName() { return std::string("CFDProjection_FreeSurface"); }

protected : 
	void buildDivergence(Grid2d&, const Grid2r& phi, const CFDsolver2D&);	
	void buildDivergence(Grid3d&, const Grid3r& phi, const CFDsolver3D&);
	
	void buildPoissonMatrix(Grid2_PoissonSolver<double>&, const Grid2r& phi, const CFDsolver2D&);
	void buildPoissonMatrix(Grid3_PoissonSolver<double>&, const Grid3r& phi, const CFDsolver3D&);	

	void buildPoissonMatrix(Grid2d& u_part, Grid2d& v_part, Grid2d& s_part, const Grid2r& phi, const CFDsolver2D&);
	void buildPoissonMatrix(Grid3d& u_part, Grid3d& v_part, Grid3d& w_part, Grid3d& s_part, const Grid3r& phi, const CFDsolver3D&);

//	void buildPoissonMatrix(Grid2_PoissonMatrix<double>&, const Grid2r& phi, const CFDsolver2D&);			
//	void buildPoissonMatrix(Grid3_PoissonMatrix<double>&, const Grid3r& phi, const CFDsolver3D&);				

	void addPressureGradient(CFDsolver2D& cfd, const Grid2d& pressure, const Grid2r& phi) { CFDProjection::addPressureGradient(cfd, pressure); }
	void addPressureGradient(CFDsolver3D& cfd, const Grid3d& pressure, const Grid3r& phi) { CFDProjection::addPressureGradient(cfd, pressure); }		
};

class NATUREDLL CFDProjection_FreeSurface_Gibou : public CFDProjection_FreeSurface {
public : 
	CFDProjection_FreeSurface_Gibou(bool tension = false, real _sigma = real(0.0728)) 
	: applySurfaceTension(tension), sigma(_sigma), maxInvTheta(20), minInvTheta(1), div2D(NULL), div3D(NULL) {}

	size_t project(CFDsolver2D& cfd, real dt, size_t iter) { return project_Liquid(cfd, *sol2D, dt, iter, thres); }
	size_t project(CFDsolver3D& cfd, real dt, size_t iter) { return project_Liquid(cfd, *sol3D, dt, iter, thres); }

	void print_log(Log& log_) const { 
		log_ << XmlElementName().c_str(); if(applySurfaceTension) log_ << " _ with Ghost Fluid Tension : " << sigma << endl; else log_ << endl;
		CFDProjection_Liquid::print_log(log_); 
	}
	static const std::string XmlElementName() { return std::string("CFDProjection_FreeSurface_Gibou"); }

protected : 
	void buildDivergence(Grid2d&, const Grid2r& phi, const CFDsolver2D&);	
	void buildDivergence(Grid3d&, const Grid3r& phi, const CFDsolver3D&);

	void buildPoissonMatrix(Grid2_PoissonSolver<double>&, const Grid2r& phi, const CFDsolver2D&);
	void buildPoissonMatrix(Grid3_PoissonSolver<double>&, const Grid3r& phi, const CFDsolver3D&);		

	void buildPoissonMatrix(Grid2d& u_part, Grid2d& v_part, Grid2d& s_part, const Grid2r& phi, const CFDsolver2D&);
	void buildPoissonMatrix(Grid3d& u_part, Grid3d& v_part, Grid3d& w_part, Grid3d& s_part, const Grid3r& phi, const CFDsolver3D&);

	void addPressureGradient(CFDsolver2D&, const Grid2d& pressure, const Grid2r& phi);
	void addPressureGradient(CFDsolver3D&, const Grid3d& pressure, const Grid3r& phi);

	const bool		applySurfaceTension;
	const double	sigma;
	const double	maxInvTheta, minInvTheta;

	Grid2d*			div2D;
	Grid3d*			div3D;

	Grid2r			curvature2D;
	Grid3r			curvature3D;
};

class NATUREDLL CFDProjection_CSF : public CFDProjection_Liquid {
public : 
	CFDProjection_CSF(double _rho_neg = 75.0f, double _rho_pos = 1.0f) : rho_neg(_rho_neg), rho_pos(_rho_pos) {}

	void set_initial_guess(Grid2_PoissonSolver<double>&) {}
	void set_initial_guess(Grid3_PoissonSolver<double>&) {}

	size_t project(CFDsolver2D& cfd, real dt, size_t iter) { return project_Liquid(cfd, *sol2D, dt, iter, thres); }
	size_t project(CFDsolver3D& cfd, real dt, size_t iter) { return project_Liquid(cfd, *sol3D, dt, iter, thres); }

	void print_log(Log& log_) const { log_ << XmlElementName().c_str() << " _ rho_neg : " << rho_neg << " _ rho_pos : " << rho_pos << endl; CFDProjection_Liquid::print_log(log_); }
	static const std::string XmlElementName() { return std::string("CFDProjection_CSF"); }

protected : 
	void buildDivergence(Grid2d& div, const Grid2r& phi, const CFDsolver2D& cfd) { CFDProjection::buildDivergence(div, cfd); }
	void buildDivergence(Grid3d& div, const Grid3r& phi, const CFDsolver3D& cfd) { CFDProjection::buildDivergence(div, cfd); }

	void buildPoissonMatrix(Grid2d& u_part, Grid2d& v_part, Grid2d& s_part, const Grid2r& phi, const CFDsolver2D&) { throwError("Not Implemented _ CFDProjection_CSF"); }
	void buildPoissonMatrix(Grid3d& u_part, Grid3d& v_part, Grid3d& w_part, Grid3d& s_part, const Grid3r& phi, const CFDsolver3D&) { throwError("Not Implemented _ CFDProjection_CSF"); }

	void buildPoissonMatrix(Grid2d& u_part, Grid2d& v_part, const Grid2r& phi, const CFDsolver2D&);
	void buildPoissonMatrix(Grid3d& u_part, Grid3d& v_part, Grid3d& w_part, const Grid3r& phi, const CFDsolver3D&);	

	void buildPoissonMatrix(Grid2_PoissonSolver<double>&, const Grid2r& phi, const CFDsolver2D&);
	void buildPoissonMatrix(Grid3_PoissonSolver<double>&, const Grid3r& phi, const CFDsolver3D&);		

	void addPressureGradient(CFDsolver2D& cfd, const Grid2d& pressure, const Grid2r& phi);
	void addPressureGradient(CFDsolver3D& cfd, const Grid3d& pressure, const Grid3r& phi);

	const double rho_neg, rho_pos;
};

};		// fluidspace

#include "CFDProjection_Project.h"

#endif