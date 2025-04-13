#include <stdafx.h>
#include <Nature/Fluid/FluidSimulator/Liquid3D.h>

using namespace fluidspace;

void Liquid3D::glBind() const {
	FluidSimulator::glBind();

	OpenGLState state = openGLState;
	state.glBind();

	glDraw();	
}

void Liquid3D::glDraw()	const {

	Material mat;
//	mat.setAmbient(RGBA::blue());
//	mat.setDiffuse(RGBA::blue());

	mat.setAmbient(RGBA::white());
	mat.setDiffuse(RGBA::white());

	if(dSurface) {
		if(dMesh)	tracking->glDrawInterface(RGBA::white(), false);
		else		tracking->glDrawInterface(RGBA::white(), true);
	}
	if(dInnerBoundary)	tracking->glDrawInnerBoundary();
	if(dBoundingBox)	solver.glDrawBoundingBox();

}

