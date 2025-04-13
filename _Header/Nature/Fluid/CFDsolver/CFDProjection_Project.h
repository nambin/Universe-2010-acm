#pragma once

#include "../../../Basis/Timer/Timer.h"

namespace fluidspace {

template<typename CFD, typename SOLVER> 
size_t CFDProjection_SinglePhase::project_SinglePhase(CFD& cfd, SOLVER& sol, size_t _iter, double _thres) {
	basis::Timer timer;	

	buildDivergence(sol.b, cfd);		
	buildPoissonMatrix(sol, cfd);		log_info << "Divergence Sum before Projection : " << sol.b.sum() << endl;

	set_initial_guess(sol);
	size_t iter = sol.solve(_iter, _thres);

	sol.x.setNeumannBoundary();
	addPressureGradient(cfd, sol.x);		

	cfd.printMaxDiv(log_info);
	log_info << "Iteration : " << iter << endl;
	timer.logTime("Single Phase Projection Time", log_info);

	if(cfd.getDimension() == 2) {	
		cfd.pressure.set(sol.x);
	}

	return iter;	
}

template<typename CFD, typename SOLVER> 
size_t CFDProjection_Liquid::project_Liquid(CFD& cfd, SOLVER& sol, real _dt, size_t _iter, double _thres) {
	dt = _dt;

	basis::Timer timer;

	if(tracking->getDimension() != cfd.getDimension())
		basis::throwError("Error in CFDProjection_FreeSurface::project()");

	typename CFD::grid_type phi; 
	tracking->getPhi(phi);

	if(liquid_biased_filter) {
		typename CFD::grid_type::value_type dx_half = phi.dx * 0.5;
		for(index_type i=0;i<phi.size();++i)
			phi[i] -= (dx_half);
	}

	cfd.printMaxDiv(sol.poissonMatrix.isBoundary, false, log_info);
	
	buildDivergence(sol.b, phi, cfd);		
	buildPoissonMatrix(sol, phi, cfd);		log_info << "Divergence Sum before Projection : " << sol.b.sum() << endl;
	
	set_initial_guess(sol);
	size_t iter = sol.solve(_iter, _thres);

	sol.x.setNeumannBoundary();	
	addPressureGradient(cfd, sol.x, phi);

	cfd.printMaxDiv(sol.poissonMatrix.isBoundary, false, log_info);
	log_info << "Iteration : " << iter << endl;
	timer.logTime("Liquid Projection Time", log_info);

	if(cfd.getDimension() == 2) {	
		cfd.pressure.set(sol.x);
	}

	return iter;
}

}