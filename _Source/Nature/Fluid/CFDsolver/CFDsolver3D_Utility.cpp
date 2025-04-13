#include <stdafx.h>
#include <Nature/Fluid/CFDsolver/CFDsolver3D.h>

using namespace std;
using namespace fluidspace;
using namespace basis;

const real CFDsolver3D::getCFL(real dt, real _h) const {	
	const MacGrid3D_VelocityGetter<real> field(u,v,w,TYPE_S);

	real max_absVel = 0;

#ifdef USE_OPENMP
	int my_end = int(nX+1);
	int i;

#pragma omp parallel for
	for( i = 1; i < my_end; ++i ) { for(index_type j=1;j<nY+1;++j) {  for(index_type k=1;k<nZ+1;++k) {

		Vector3<real> vel = field.getVel(i,j,k);
		real absVel = abs(vel[0]) + abs(vel[1]) + abs(vel[2]);

		if(max_absVel < absVel) {
#pragma omp critical
			{
				max_absVel = (absVel > max_absVel) ? absVel : max_absVel;
			}
		}
	}}}

#else
	for(index_type i=1;i<nX+1;++i) { for(index_type j=1;j<nY+1;++j) {  for(index_type k=1;k<nZ+1;++k) {

		Vector3<real> vel = field.getVel(i,j,k);
		real absVel = abs(vel[0]) + abs(vel[1]) + abs(vel[2]);

		if(max_absVel < absVel) 
			max_absVel = absVel;
	} } }
#endif

	return dt * max_absVel / _h;
}
