#include <stdafx.h>
#include <Nature/Fluid/FluidSimulator/Smoke3D.h>
#include <Principle/VirtualSpace/Material/Material.h>
#include <Principle/VirtualSpace/Light/Light.h>

using namespace fluidspace;

void Smoke3D::glBind() const {
	FluidSimulator::glBind();

	OpenGLState state = openGLState;
	state.glBind();

	glDraw();	
}

void Smoke3D::glDraw() const {
	if(dBoundingBox)	densityField.enableDrawBoundingBox();
	else				densityField.disableDrawBoundingBox();

	if(dDensity)		densityField.glDraw();
}

void Smoke3D::resizeScreen(unsigned int width, unsigned int height) {
	densityField.initialize_Shader(width, height);
}