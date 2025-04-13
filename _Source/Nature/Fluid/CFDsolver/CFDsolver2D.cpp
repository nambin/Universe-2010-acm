#include <stdafx.h>
#include <Nature/Fluid/CFDsolver/CFDsolver2D.h>
#include <Nature/Fluid/FluidScene/FluidScene_Objects.h>

using namespace std;
using namespace fluidspace;

void CFDsolver2D::set(index_type _nX, index_type _nY, real _h)	{
	
	MacGrid2D::set(_nX, _nY, _h);

	const Vector2<real> cOffset(-0.5f, -0.5f), uOffset(0.0f, -0.5f), vOffset(-0.5f, 0.0f);

	isKinematicBoundary.init(nX+2, nY+2, h, h, cOffset);
	isKinematicBoundary.zero();
	isKinematicBoundary.set_boundary(true);

	u.init(nX+1, nY+2, h, h, uOffset);	u.zero();
	u0.init(nX+1, nY+2, h, h, uOffset);	u0.zero();
	v.init(nX+2, nY+1, h, h, vOffset);	v.zero();
	v0.init(nX+2, nY+1, h, h, vOffset);	v0.zero();

	pressure.init(nX+2, nY+2, h, h, cOffset);	pressure.zero();

	log_info << "CFDSolver2D::set()" << endl;
	log_info << "(nX, nY, h) = " << nX << ", " << nY << ", " << h << endl;
	log_info << "(Lx, Ly)    = " << nX*h << ", " << nY*h << endl;
}

void CFDsolver2D::setVelocityField(const VelocityField2D<real> &field) {
	for( index_type i=0; i<u.nX; ++i ) 
		for( index_type j=0; j<u.nY; ++j)	
			u(i,j) = u0(i,j) = field.getU(u.pos(i,j));
	
	for( index_type i=0; i<v.nX; ++i) 
		for( index_type j=0; j<v.nY; ++j)	
			v(i,j) = v0(i,j) = field.getV(v.pos(i,j));
}

void CFDsolver2D::update_InternalBoundary(const FluidScene_Objects& _objs) {

	if(_objs.getDimension() != 2)
		throwError("Error in CFDsolver2D::update_InternalBoundary(const FluidScene_Objects& objs)");

	const FluidScene_Objects2D& objs = static_cast<const FluidScene_Objects2D&>(_objs);

	const real dx_half = isKinematicBoundary.dx * 0.5f;
	const real dy_half = isKinematicBoundary.dy * 0.5f;

	isKinematicBoundary.set_boundary(true);

	for(index_type i=1;i<isKinematicBoundary.nX-1;++i) {
		for(index_type j=1;j<isKinematicBoundary.nY-1;++j) {
			isKinematicBoundary(i,j) = objs.isInside(isKinematicBoundary.pos(i,j)); 

			//isKinematicBoundary(i,j) =	(  objs.isInside(isKinematicBoundary.pos(i,j) + Vector2<real>( dx_half, dy_half)) 
			//							|| objs.isInside(isKinematicBoundary.pos(i,j) + Vector2<real>(-dx_half, dy_half)) 
			//							|| objs.isInside(isKinematicBoundary.pos(i,j) + Vector2<real>( dx_half,-dy_half)) 
			//							|| objs.isInside(isKinematicBoundary.pos(i,j) + Vector2<real>(-dx_half,-dy_half)) );
		}
	}
}

void CFDsolver2D::interp_KinematicBoundary(Grid2b& isbdry, const MacGrid2D<real>::FIELD_TYPE type) const {
	switch(type) {
		case TYPE_S :	isbdry = isKinematicBoundary;
						break;

		case TYPE_U :	isbdry.set_boundary(true);
						for(index_type i=1;i<isbdry.nX-1;++i) 
							for(index_type j=1;j<isbdry.nY-1;++j) 
								isbdry(i,j) = (isKinematicBoundary(i,j) || isKinematicBoundary(i+1,j));
						break;

		case TYPE_V :	isbdry.set_boundary(true);
						for(index_type i=1;i<isbdry.nX-1;++i) 
							for(index_type j=1;j<isbdry.nY-1;++j) 
								isbdry(i,j) = (isKinematicBoundary(i,j) || isKinematicBoundary(i,j+1));
						break;

		default : throwError("Not Implemented");
	}	
}

void CFDsolver2D::setVelocity_InternalBoundary(const FluidScene_Objects* _objects) {
	if(_objects == NULL)
		return;

	if(_objects->getDimension() != 2)
		throwError("Error in CFDsolver2D::setVelocity_InternalBoundary(const FluidScene_Objects* objects)");

	const FluidScene_Objects2D* objects = static_cast<const FluidScene_Objects2D*>(_objects);

	for(index_type i=1;i<u.nX-1;++i) {			
		for(index_type j=1;j<u.nY-1;++j) {
			if(isKinematicBoundary(i,j) || isKinematicBoundary(i+1,j)) {
				u(i,j) = objects->getVel(u.pos(i,j))[0];
			}
		}}

	for(index_type i=1;i<v.nX-1;++i) {		
		for(index_type j=1;j<v.nY-1;++j) {
			if(isKinematicBoundary(i,j) || isKinematicBoundary(i,j+1)) {
				v(i,j) = objects->getVel(v.pos(i,j))[1];
			}
		}}
}