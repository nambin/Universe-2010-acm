#include <stdafx.h>
#include <Nature/Fluid/FluidSimulator/Liquid2D.h>

using namespace fluidspace;

void Liquid2D::glBind() const {
	FluidSimulator::glBind();

	OpenGLState state = openGLState;
	state.disable_Lighting();
	state.disable_ColorMaterial();
	state.glBind();

	glDraw();	
}

void Liquid2D::glDraw() const {

	if(dLiquid)	{
		tracking->glDrawInterface(RGBA::blue(), dFilledInterface);
	}
	if(dInternalBoundary) {
		const InternalBoundary* internalBdry = integrator->integrator->internalBdry;
		if(internalBdry)
			internalBdry->glDraw();
	}
	if(dPhi)			
		tracking->glDrawPhi();
	if(dAuxiliary) {
		InterfaceTracking_PLS2D* pls = dynamic_cast<InterfaceTracking_PLS2D*>(tracking);
		if(pls != NULL)	pls->glDrawParticles();		

		InterfaceTracking_SRL2D* srl = dynamic_cast<InterfaceTracking_SRL2D*>(tracking);
		if(srl != NULL && dPhi == false)	srl->glDrawAuxiliary();
		if(srl != NULL && dPhi == true)		srl->glDrawPhi_SubGrid();
	}
	if(dVelocity)		solver.glDrawVelocity(scaleVelocity);
	if(dPressure)		solver.glDrawPressure();
	if(dGrid)			solver.glDrawUniformGrid(RGBA::white() * 0.3f, 0.5f);
	if(dBoundingBox)	solver.glDrawBoundingBox(RGBA::white(), 3.0f);

}