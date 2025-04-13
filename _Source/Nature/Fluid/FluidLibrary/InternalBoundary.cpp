#include <stdafx.h>
#include <Nature/Fluid/FluidLibrary/InternalBoundary.h>
#include <Nature/Fluid/InterfaceTracking/SRL_Extrapolation.h>
#include <Principle/Geometry/Levelset/Levelset_Redistancing.h>

using namespace std;
using namespace fluidspace;

void InternalBoundary_Implicit2D::update_InternalBoundary(const FluidScene_Objects* _objs, const InterfaceTracking& tracking, const CFDsolver& _cfd) {

	const FluidScene_Objects2D* objs = static_cast<const FluidScene_Objects2D*>(_objs);
	const CFDsolver2D& cfd = static_cast<const CFDsolver2D&>(_cfd);

	phi.init(tracking.get_nX(), tracking.get_nY(), tracking.get_h(), tracking.get_h(), tracking.get_offset2());	

	const real half_h = phi.dx * 0.5f;

	if(phi.nX == cfd.isKinematicBoundary.nX && phi.nY == cfd.isKinematicBoundary.nY) {
		for( grid_index_type i = 0; i < phi.nX; ++i ) 
			for(grid_index_type j=0;j<phi.nY;++j)
				phi(i,j) = (cfd.isKinematicBoundary(i,j)) ? half_h : -half_h;
	}
	else {
		for( grid_index_type i = 0; i < phi.nX; ++i ) 
			for(grid_index_type j=0;j<phi.nY;++j)
				phi(i,j) = (cfd.isKinematicBoundary.get_nearest(phi.pos(i,j))) ? half_h : -half_h;
	}

	Levelset_Redistancing_PDE_1stUpwind<real> redist(band, 0.5);
	redist.redistance(phi);
}

void InternalBoundary_Implicit2D::extrapolate_phi(InterfaceTracking& tracking, real band) {
	if(dynamic_cast<InterfaceTracking_Levelset2D*>(&tracking) != NULL) extrapolate_phi(static_cast<InterfaceTracking_Levelset2D&>(tracking), band);
	else if(dynamic_cast<InterfaceTracking_SRL2D*>(&tracking) != NULL) extrapolate_phi(static_cast<InterfaceTracking_SRL2D&>(tracking), band);
}

void InternalBoundary_Implicit2D::extrapolate_phi(InterfaceTracking_Levelset2D& lv, real band) {
	Levelset_Extrapolation_PDE<real> extrapolation(band);
	extrapolation.extrapolate(lv.phi, MacGrid2D<real>::TYPE_S, GridBoundaryCondition_Neumann<real,real>(), phi);
}

void InternalBoundary_Implicit2D::extrapolate_phi(InterfaceTracking_SRL2D& srl, real band) {
	Levelset_Extrapolation_PDE<real> extrapolation(band);
	extrapolation.extrapolate(srl.srl->phi, MacGrid2D<real>::TYPE_S, GridBoundaryCondition_Neumann<real,real>(), phi);

	SRL_Extrapolation_PDE srl_extrapolation;
	srl_extrapolation.extrapolate(srl, phi, 0);
	srl.srl->synchronize_grid_values();
}

void InternalBoundary_Implicit2D::glDraw() const {
	const RGBA cp = RGBA::brown();
	const RGBA cn = RGBA::white();
	const real band = 7.0f, dx_half = phi.dx*0.5f, dy_half = phi.dy*0.5f;	

	real eps = 0.003f;

	Vector3f p0, p1, p2, p3;
	for(grid_index_type i=1;i<phi.nX-1;i++)	{
		for(grid_index_type j=1;j<phi.nY-1;j++)	{
			
			if( phi(i,j) < 0 ) continue;

			const real scale = (band - abs(phi(i,j)) * phi.dxInv) / band;
			const Vector2<real> pos = phi.pos(i,j);

			if(scale < 0)			continue;
			else if(phi(i,j) > 0)	glColor4f(cp[0] * scale, cp[1] * scale, cp[2] * scale, cp[3]);
			else					glColor4f(cn[0] * scale, cn[1] * scale, cn[2] * scale, cn[3]);

			p0 = Vector3f(pos[0]-dx_half, pos[1]-dy_half, eps);
			p1 = Vector3f(pos[0]+dx_half, pos[1]-dy_half, eps);
			p2 = Vector3f(pos[0]+dx_half, pos[1]+dy_half, eps);
			p3 = Vector3f(pos[0]-dx_half, pos[1]+dy_half, eps);			

			glLineWidth(3.0f);
			if( phi(i+1, j) < 0) { glBegin(GL_LINES); glVertex3f(p1[0], p1[1], p1[2]); glVertex3f(p2[0], p2[1], p2[2]); glEnd(); }
			if( phi(i-1, j) < 0) { glBegin(GL_LINES); glVertex3f(p0[0], p0[1], p0[2]); glVertex3f(p3[0], p3[1], p3[2]); glEnd(); }
			if( phi(i, j+1) < 0) { glBegin(GL_LINES); glVertex3f(p2[0], p2[1], p2[2]); glVertex3f(p3[0], p3[1], p3[2]); glEnd(); }
			if( phi(i, j-1) < 0) { glBegin(GL_LINES); glVertex3f(p1[0], p1[1], p1[2]); glVertex3f(p0[0], p0[1], p0[2]); glEnd(); }
			
			//glBegin(GL_QUADS);
			//	glVertex3f(p0[0], p0[1], p0[2]);
			//	glVertex3f(p1[0], p1[1], p1[2]);
			//	glVertex3f(p2[0], p2[1], p2[2]);
			//	glVertex3f(p3[0], p3[1], p3[2]);
			//glEnd();			
		}
	}	
}

void InternalBoundary_Implicit3D::update_InternalBoundary(const FluidScene_Objects* _objs, const InterfaceTracking& tracking, const CFDsolver& _cfd) {
	const FluidScene_Objects3D* objs = static_cast<const FluidScene_Objects3D*>(_objs);
	const CFDsolver3D& cfd = static_cast<const CFDsolver3D&>(_cfd);

	phi.init(tracking.get_nX(), tracking.get_nY(), tracking.get_nZ(), tracking.get_h(), tracking.get_h(), tracking.get_h(), tracking.get_offset3());	

	const real half_h = phi.dx * 0.5f;

	if(phi.nX == cfd.isKinematicBoundary.nX && phi.nY == cfd.isKinematicBoundary.nY && phi.nZ == cfd.isKinematicBoundary.nZ) {

#ifdef USE_OPENMP
		int my_end = int(phi.nX);
		int i;

#pragma omp parallel for
#else
		grid_index_type my_end = phi.nX;
		grid_index_type i;
#endif
		for( i = 0; i < my_end; ++i ) 
			for(grid_index_type j=0;j<phi.nY;++j)
				for(grid_index_type k=0;k<phi.nZ;++k)
					phi(i,j,k) = (cfd.isKinematicBoundary(i,j,k)) ? half_h : -half_h;
	}
	else {

#ifdef USE_OPENMP
		int my_end = int(phi.nX);
		int i;

#pragma omp parallel for
#else
		grid_index_type my_end = phi.nX;
		grid_index_type i;
#endif
		for( i = 0; i < my_end; ++i ) 
			for(grid_index_type j=0;j<phi.nY;++j)
				for(grid_index_type k=0;k<phi.nZ;++k)
					phi(i,j,k) = (cfd.isKinematicBoundary.get_nearest(phi.pos(i,j,k))) ? half_h : -half_h;
	}

	Levelset_Redistancing_PDE_1stUpwind<real> redist(band, 0.5);
	redist.redistance(phi);
}

void InternalBoundary_Implicit3D::extrapolate_phi(InterfaceTracking& tracking, real band) {
	if(dynamic_cast<InterfaceTracking_Levelset3D*>(&tracking) != NULL) extrapolate_phi(static_cast<InterfaceTracking_Levelset3D&>(tracking), band);
	else if(dynamic_cast<InterfaceTracking_SRL3D*>(&tracking) != NULL) extrapolate_phi(static_cast<InterfaceTracking_SRL3D&>(tracking), band);
}

void InternalBoundary_Implicit3D::extrapolate_phi(InterfaceTracking_Levelset3D& lv, real band) {
	Levelset_Extrapolation_PDE<real> extrapolation(band);
	extrapolation.extrapolate(lv.phi, MacGrid3D<real>::TYPE_S, GridBoundaryCondition_Neumann<real,real>(), phi);
}

void InternalBoundary_Implicit3D::extrapolate_phi(InterfaceTracking_SRL3D& srl, real band) {
	Levelset_Extrapolation_PDE<real> extrapolation(band);
	extrapolation.extrapolate(srl.srl->phi, MacGrid3D<real>::TYPE_S, GridBoundaryCondition_Neumann<real,real>(), phi);

	SRL_Extrapolation_PDE srl_extrapolation;
	srl_extrapolation.extrapolate(srl, phi, 0);
	srl.srl->synchronize_grid_values();
}

void InternalBoundary_Implicit3D::glDraw() const {}
