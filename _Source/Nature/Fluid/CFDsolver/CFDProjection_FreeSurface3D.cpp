#include <stdafx.h>
#include <Nature/Fluid/CFDsolver/CFDProjection.h>

using namespace std;
using namespace fluidspace;

//size_t CFDProjection_FreeSurface::project(CFDsolver3D& cfd, size_t _iter, real _thres) {
//	Timer timer;
//
//	if(tracking->getDimension() != 3)
//		basis::throwError("Error in CFDProjection_FreeSurface::project()");
//
//	Grid3r phi3D;
//	tracking->getPhi(phi3D);
//
//	buildDivergenceFreeSurface(sol3D->b, phi3D, cfd);
//	buildPoissonMatrix(*sol3D, phi3D, cfd);	
//
//	size_t iter = sol3D->solve(_iter, _thres);
//
//	sol3D->x.setNeumannBoundary();
//	CFDProjection::addPressureGradient(cfd, sol3D->x);
//
//	return iter;
//}
//
