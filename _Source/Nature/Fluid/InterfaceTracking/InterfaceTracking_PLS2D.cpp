#include <stdafx.h>
#include <Nature/Fluid/InterfaceTracking/InterfaceTracking_PLS2D.h>

using namespace std;
using namespace fluidspace;

void InterfaceTracking_PLS2D::set(index_type nX, index_type nY, real h) {
	InterfaceTracking_Levelset2D::set(nX, nY, h);

	log_info << "InterfaceTracking_PLS2D::set()" << endl;
	log_info << "(nX, nY, h) = " << nX << ", " << nY << ", " << h << endl;
}

void InterfaceTracking_PLS2D::setPhi(const Grid2r &_phi) {
	InterfaceTracking_Levelset2D::setPhi(_phi);
}

void InterfaceTracking_PLS2D::setPhi(const geometryspace::DistanceField2D<real> &field) {
	InterfaceTracking_Levelset2D::setPhi(field);
}

void InterfaceTracking_PLS2D::glDrawParticles(const RGBA &color0, const RGBA &color1, GLfloat pointsize) const {

	glPointSize(pointsize);

	glBegin(GL_POINTS);
	for(size_t i=0;i<particles.size();++i) {
		if(particles[i].isDeleted)
			continue;

		if(particles[i].sign == 1)	glColor4f(color0[0], color0[1], color0[2], color0[3]);
		else						glColor4f(color1[0], color1[1], color1[2], color1[3]);

		glVertex3f(particles[i].pos[0], particles[i].pos[1], 0.0005f);
	}	
	glEnd();
}
