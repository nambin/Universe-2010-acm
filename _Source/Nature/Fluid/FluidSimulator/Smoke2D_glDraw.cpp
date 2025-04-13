#include <stdafx.h>
#include <Nature/Fluid/FluidSimulator/Smoke2D.h>
#include <Principle/VirtualSpace/Material/Material.h>
#include <Principle/VirtualSpace/Light/Light.h>

using namespace fluidspace;

void Smoke2D::glBind() const {
	FluidSimulator::glBind();

	OpenGLState state = openGLState;
	state.glBind();

	glDisable(GL_COLOR_MATERIAL);
	glDisable(GL_LIGHTING);

	glDraw();	
}

void Smoke2D::glDraw() const {

	if(dKinematicBoundary)	solver.glDrawKinematicBoundary();
	if(dPressure)			solver.glDrawPressure();
	if(dVelocity)			solver.glDrawVelocity(scaleVelocity);
	if(dBoundingBox)		solver.glDrawBoundingBox(RGBA::white(), 3.0f);
	if(dDensity)			densityField.glDraw(RGBA::white());		

}

