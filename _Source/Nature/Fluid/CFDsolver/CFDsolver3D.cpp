#include <stdafx.h>
#include <Nature/Fluid/CFDsolver/CFDsolver3D.h>
#include <Nature/Fluid/FluidScene/FluidScene_Objects.h>

using namespace std;
using namespace fluidspace;

void CFDsolver3D::set(index_type _nX, index_type _nY, index_type _nZ, real _h)	{
	MacGrid3D::set(_nX, _nY, _nZ, _h);

	const Vector3<real> cOffset(-0.5f, -0.5f, -0.5f), uOffset(0.0f, -0.5f, -0.5f), vOffset(-0.5f, 0.0f, -0.5f), wOffset(-0.5f, -0.5f, 0.0f);

	isKinematicBoundary.init(nX+2, nY+2, nZ+2, h, h, h, cOffset);	
	isKinematicBoundary.zero();
	isKinematicBoundary.set_boundary(true);
	
	u.init(nX+1, nY+2, nZ+2, h, h, h, uOffset);		u.zero();
	u0.init(nX+1, nY+2, nZ+2, h, h, h, uOffset);	u0.zero();
	v.init(nX+2, nY+1, nZ+2, h, h, h, vOffset);		v.zero();
	v0.init(nX+2, nY+1, nZ+2, h, h, h, vOffset);	v0.zero();
	w.init(nX+2, nY+2, nZ+1, h, h, h, wOffset);		w.zero();
	w0.init(nX+2, nY+2, nZ+1, h, h, h, wOffset);	w0.zero();

	log_info << "CFDSolver3D::set()" << endl;
	log_info << "(nX, nY, nZ, h) = " << nX << ", " << nY << ", " << nZ << ", " << h << endl;
	log_info << "(Lx, Ly, Lz)    = " << nX*h << ", " << nY*h << ", " << nZ*h << endl;
}

void CFDsolver3D::setVelocityField(const VelocityField3D<real> &field) {
	for( index_type i=0; i<u.nX; ++i ) 
		for( index_type j=0; j<u.nY; ++j)	
			for( index_type k=0; k<u.nZ; ++k)	
				u(i,j,k) = u0(i,j,k) = field.getU(u.pos(i,j,k));

	for( index_type i=0; i<v.nX; ++i ) 
		for( index_type j=0; j<v.nY; ++j)	
			for( index_type k=0; k<v.nZ; ++k)	
				v(i,j,k) = v0(i,j,k) = field.getV(v.pos(i,j,k));

	for( index_type i=0; i<w.nX; ++i ) 
		for( index_type j=0; j<w.nY; ++j)	
			for( index_type k=0; k<w.nZ; ++k)	
				w(i,j,k) = w0(i,j,k) = field.getW(w.pos(i,j,k));
}

void CFDsolver3D::update_InternalBoundary(const FluidScene_Objects& _objs) {

	if(_objs.getDimension() != 3)
		throwError("Error in CFDsolver2D::update_InternalBoundary(const FluidScene_Objects& objs)");

	const FluidScene_Objects3D& objs = static_cast<const FluidScene_Objects3D&>(_objs);

	isKinematicBoundary.set_boundary(true);

#ifdef USE_OPENMP
	int i;
	int my_end = int(isKinematicBoundary.nX-1);	
#pragma omp parallel for
#else
	index_type i;
	index_type my_end = isKinematicBoundary.nX-1;	
#endif
	for( i = 1; i < my_end; ++i ) {	
		for(index_type j=1;j<isKinematicBoundary.nY-1;++j) {
			for(index_type k=1;k<isKinematicBoundary.nZ-1;++k) {
				isKinematicBoundary(i,j,k) = objs.isInside(isKinematicBoundary.pos(i,j,k));
			}
		}
	}
}

void CFDsolver3D::interp_KinematicBoundary(Grid3b& isbdry, const MacGrid3D<real>::FIELD_TYPE type) const {
	switch(type) {
		case TYPE_S :	isbdry = isKinematicBoundary;
						break;

		case TYPE_U : {
						isbdry.set_boundary(true);

				#ifdef USE_OPENMP
						int my_end = int(isbdry.nX-1);
						int i;

				#pragma omp parallel for
				#else
						index_type my_end = isbdry.nX-1;
						index_type i;
				#endif
						for( i = 1; i < my_end; ++i ) 					
							for(index_type j=1;j<isbdry.nY-1;++j) 
								for(index_type k=1;k<isbdry.nZ-1;++k) 
									isbdry(i,j,k) = (isKinematicBoundary(i,j,k) || isKinematicBoundary(i+1,j,k));
						break;
					  }

		case TYPE_V : {
						isbdry.set_boundary(true);

				#ifdef USE_OPENMP
						int my_end = int(isbdry.nX-1);
						int i;

				#pragma omp parallel for
				#else
						index_type my_end = isbdry.nX-1;
						index_type i;
				#endif
						for( i = 1; i < my_end; ++i ) 					
							for(index_type j=1;j<isbdry.nY-1;++j) 
								for(index_type k=1;k<isbdry.nZ-1;++k) 
									isbdry(i,j,k) = (isKinematicBoundary(i,j,k) || isKinematicBoundary(i,j+1,k));
						break;
					  }

		case TYPE_W : {
						isbdry.set_boundary(true);

				#ifdef USE_OPENMP
						int my_end = int(isbdry.nX-1);
						int i;

				#pragma omp parallel for
				#else
						index_type my_end = isbdry.nX-1;
						index_type i;
				#endif
						for( i = 1; i < my_end; ++i ) 					
							for(index_type j=1;j<isbdry.nY-1;++j) 
								for(index_type k=1;k<isbdry.nZ-1;++k) 
									isbdry(i,j,k) = (isKinematicBoundary(i,j,k) || isKinematicBoundary(i,j,k+1));
						break;
					  }

		default : throwError("Not Implemented");
	}	
}

void CFDsolver3D::setVelocity_InternalBoundary(const FluidScene_Objects* _objects) {
	if(_objects == NULL)
		return;

	if(_objects->getDimension() != 3)
		throwError("Error in CFDsolver3D::setVelocity_InternalBoundary(const FluidScene_Objects* objects)");

	const FluidScene_Objects3D* objects = static_cast<const FluidScene_Objects3D*>(_objects);

#ifdef USE_OPENMP
	int my_end = int(u.nX-1);
	int i;

#pragma omp parallel for
#else
	index_type my_end = u.nX-1;
	index_type i;
#endif
	for( i = 1; i < my_end; ++i ) {		
		for(index_type j=1;j<u.nY-1;++j) {
			for(index_type k=1;k<u.nZ-1;++k) {
				if(isKinematicBoundary(i,j,k) || isKinematicBoundary(i+1,j,k)) {
					u(i,j,k) = objects->getVel(u.pos(i,j,k))[0];
				}
			}}}

#ifdef USE_OPENMP
	my_end = int(v.nX-1);
#pragma omp parallel for
#else
	my_end = v.nX-1;
#endif
	for( i = 1; i < my_end; ++i ) {		
		for(index_type j=1;j<v.nY-1;++j) {
			for(index_type k=1;k<v.nZ-1;++k) {
				if(isKinematicBoundary(i,j,k) || isKinematicBoundary(i,j+1,k)) {
					v(i,j,k) = objects->getVel(v.pos(i,j,k))[1];
				}
			}}}

#ifdef USE_OPENMP
	my_end = int(w.nX-1);
#pragma omp parallel for
#else
	my_end = w.nX-1;
#endif
	for( i = 1; i < my_end; ++i ) {		
		for(index_type j=1;j<w.nY-1;++j) {
			for(index_type k=1;k<w.nZ-1;++k) {
				if(isKinematicBoundary(i,j,k) || isKinematicBoundary(i,j,k+1)) {
					w(i,j,k) = objects->getVel(w.pos(i,j,k))[2];
				}
			}}}
}

void CFDsolver3D::printDivergenceInfo(index_type i, index_type j, index_type k) const {
	log_info << i << " " << j << " " << k << " Divergence Info" << endl;

	log_info << "u -1 _ " << u(iS2U(i,j,k,0)) << " "  << "u +1 _ " << u(iS2U(i,j,k,1)) << endl;
	log_info << "v -1 _ " << v(iS2V(i,j,k,0)) << " "  << "v +1 _ " << v(iS2V(i,j,k,1)) << endl;
	log_info << "w -1 _ " << w(iS2W(i,j,k,0)) << " "  << "w +1 _ " << w(iS2W(i,j,k,1)) << endl;

	log_info << "Divergence : " << (  u(iS2U(i,j,k,1)) - u(iS2U(i,j,k,0)) 
									+ v(iS2V(i,j,k,1)) - v(iS2V(i,j,k,0)) 
									+ w(iS2W(i,j,k,1)) - w(iS2W(i,j,k,0))) * hInv << endl << endl;
}