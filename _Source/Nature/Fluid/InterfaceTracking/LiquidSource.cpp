#include <stdafx.h>
#include <Nature/Fluid/InterfaceTracking/LiquidSource.h>

using namespace std;
using namespace fluidspace;

void LiquidSource::setVel(CFDsolver2D& cfd, const Vector2<real>& vel, const geometryspace::DistanceField2D<real>& field, const real _eps) const {
	const real eps = _eps * cfd.h;

	for(CFDsolver2D::index_type i=0;i<cfd.u.nX;++i) 
		for(CFDsolver2D::index_type j=0;j<cfd.u.nY;++j) 
			if(field.getPhi(cfd.u.pos(i,j)) < eps)
				cfd.u(i,j) = vel[0];

	for(CFDsolver2D::index_type i=0;i<cfd.v.nX;++i) 
		for(CFDsolver2D::index_type j=0;j<cfd.v.nY;++j) 
			if(field.getPhi(cfd.v.pos(i,j)) < eps)
				cfd.v(i,j) = vel[1];
}

void LiquidSource::setVel(CFDsolver3D& cfd, const Vector3<real>& vel, const geometryspace::DistanceField3D<real>& field, const real _eps) const {
	const real eps = _eps * cfd.h;

	for(CFDsolver3D::index_type i=0;i<cfd.u.nX;++i) {
		for(CFDsolver3D::index_type j=0;j<cfd.u.nY;++j) {
			for(CFDsolver3D::index_type k=0;k<cfd.u.nZ;++k) {
				if(field.getPhi(cfd.u.pos(i,j,k)) < eps) {
					cfd.u(i,j,k) = vel[0];
				}
			}}}

	for(CFDsolver3D::index_type i=0;i<cfd.v.nX;++i) {
		for(CFDsolver3D::index_type j=0;j<cfd.v.nY;++j) {
			for(CFDsolver3D::index_type k=0;k<cfd.v.nZ;++k) {
				if(field.getPhi(cfd.v.pos(i,j,k)) < eps) {
					cfd.v(i,j,k) = vel[1];

					if( j < 10 )
						cout << "xxx";
				}
			}}}


	for(CFDsolver3D::index_type i=0;i<cfd.w.nX;++i) {
		for(CFDsolver3D::index_type j=0;j<cfd.w.nY;++j) {
			for(CFDsolver3D::index_type k=0;k<cfd.w.nZ;++k) {
				if(field.getPhi(cfd.w.pos(i,j,k)) < eps) {
					cfd.w(i,j,k) = vel[2];
				}
			}}}
}

void LiquidSource_Continuous::add(InterfaceTracking* tracking, CFDsolver& cfd, real dt) const {

	if(field->getDimension() != tracking->getDimension() || field->getDimension() != cfd.getDimension())
		throwError("Error in LiquidSource_Continuous::add()");

	for(size_t i=0;i<interval.size();i+=2) {
		if(interval[i] <= time && time < interval[i+1]) {
			if(tracking->getDimension() == 2) {
				DistanceField2D<real>& field2D = (DistanceField2D<real>&)(*field);
				tracking->addSource(field2D);

				CFDsolver2D& cfd2D = (CFDsolver2D&)cfd;
				setVel(cfd2D, vel2D, field2D);				
			}

			if(tracking->getDimension() == 3) {
				DistanceField3D<real>& field3D = (DistanceField3D<real>&)(*field);
				tracking->addSource(field3D);

				CFDsolver3D& cfd3D = (CFDsolver3D&)cfd;
				setVel(cfd3D, vel3D, field3D);			
			}
			break;
		}
	}

	time += dt;
}

void LiquidSource_Continuous::add_interval(fluidspace::real t0, fluidspace::real t1) {
	interval.push_back(t0);
	interval.push_back(t1);
}

void LiquidSource_Discrete::add(InterfaceTracking* tracking, CFDsolver& cfd, real dt) const {

	if(field->getDimension() != tracking->getDimension())
		throwError("Error in LiquidSource_Continuous::add()");

	for(size_t i=0;i<time_point.size();i++) {
		if(time <= time_point[i] && time_point[i] < time + dt) {
			if(tracking->getDimension() == 2) {
				DistanceField2D<real>& field2D = (DistanceField2D<real>&)(*field);
				tracking->addSource(field2D);

				CFDsolver2D& cfd2D = (CFDsolver2D&)cfd;
				setVel(cfd2D, vel2D, field2D);				
			}

			if(tracking->getDimension() == 3) {				
				DistanceField3D<real>& field3D = (DistanceField3D<real>&)(*field);
				tracking->addSource(field3D);

				CFDsolver3D& cfd3D = (CFDsolver3D&)cfd;
				setVel(cfd3D, vel3D, field3D);		
			}
			break;
		}
	}

	time += dt;
}

void LiquidSource_Discrete::add_time_point(real t0) {
	time_point.push_back(t0);
	sort(time_point.begin(), time_point.end());
}