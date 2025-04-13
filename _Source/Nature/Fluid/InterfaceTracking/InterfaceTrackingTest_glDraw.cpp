#include <stdafx.h>
#include <Nature/Fluid/InterfaceTracking/InterfaceTrackingTest.h>

using namespace fluidspace;

void InterfaceTrackingTest::glBind() const {
	FluidSimulator::glBind();

	OpenGLState state = openGLState;
	state.disable_Lighting();
	state.disable_ColorMaterial();
	state.glBind();

	glDraw();	
}

void InterfaceTrackingTest::glDraw() const {

	if(tracking == NULL) return;

	Material mat;
	mat.setAmbient(RGBA::blue());
	mat.setDiffuse(RGBA::blue());
	mat.setSpecular(RGBA::blue() * 0.25f);
	mat.setShininess(8.5f);

	if(dLiquid) { 
		tracking->glDrawInterface(mat, dFilledMesh);	
	} 
	if(dPhi && tracking->getDimension() == 2 && tracking->isEulerian())	{ tracking->glDrawPhi();				}

	if(dAuxiliary) {
		InterfaceTracking_PLS2D* pls = dynamic_cast<InterfaceTracking_PLS2D*>(tracking);
		if(pls != NULL)	pls->glDrawParticles();		

		InterfaceTracking_SRL2D* srl = dynamic_cast<InterfaceTracking_SRL2D*>(tracking);
		if(srl != NULL && dPhi == false)	srl->glDrawAuxiliary();
		if(srl != NULL && dPhi == true)		srl->glDrawPhi_SubGrid();
	}

	if(dGrid && tracking->getDimension() == 2 && tracking->isEulerian())	{ 
		tracking->glDrawGrid();			
	}
	if(dBoundingBox && tracking->isEulerian())								{ 
		tracking->glDrawBoundingBox();	
	}

}