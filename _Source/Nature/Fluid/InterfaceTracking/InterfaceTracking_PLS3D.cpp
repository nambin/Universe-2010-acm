#include <stdafx.h>
#include <Nature/Fluid/InterfaceTracking/InterfaceTracking_PLS3D.h>

using namespace std;
using namespace fluidspace;

void InterfaceTracking_PLS3D::set(index_type nX, index_type nY, index_type nZ, real h) {
	InterfaceTracking_Levelset3D::set(nX, nY, nZ, h);

	log_info << "InterfaceTracking_PLS3D::set()" << endl;
	log_info << "(nX, nY, nZ, h) = " << nX << ", " << nY << ", " << nZ << ", " << h << endl;
}

void InterfaceTracking_PLS3D::setPhi(const Grid3r &_phi) {
	InterfaceTracking_Levelset3D::setPhi(_phi);
}

void InterfaceTracking_PLS3D::setPhi(const geometryspace::DistanceField3D<real> &field) {
	InterfaceTracking_Levelset3D::setPhi(field);
}

void InterfaceTracking_PLS3D::glDrawParticles(const RGBA &color0, const RGBA &color1, GLfloat pointsize) const {

	glPointSize(pointsize);

	glBegin(GL_POINTS);
	for(size_t i=0;i<particles.size();++i) {
		if(particles[i].isDeleted)
			continue;

		if(particles[i].sign == 1)	glColor4f(color0[0], color0[1], color0[2], color0[3]);
		else						glColor4f(color1[0], color1[1], color1[2], color1[3]);

		glVertex3f(particles[i].pos[0], particles[i].pos[1], particles[i].pos[2]);
	}	
	glEnd();
}
