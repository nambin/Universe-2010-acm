#include <stdafx.h>
#include <Principle/VirtualSpace/OpenGLState/OpenGLState.h>

using namespace virtualspace;

OpenGLState::OpenGLState() {
	init();
}

OpenGLState::OpenGLState(const OpenGLState &state) {
	init();
	(*this) = state;
}

void OpenGLState::operator=(const OpenGLState &state) {
	for(unsigned int i=0;i<binaryState.size();++i) 
		*(binaryState[i]) = *(state.binaryState[i]);

	for(unsigned int i=0;i<twoInputState.size();++i) 
		*(twoInputState[i]) = *(state.twoInputState[i]);
}

void OpenGLState::init() {
	binaryState.resize(numBinaryState);

	OpenGLBinaryStateVariable *state = NULL;

	state = new OpenGLBinaryStateLighting();		binaryState[IDX_BINARY_LIGHTING] = state;
	state = new OpenGLBinaryStateDepthTest();		binaryState[IDX_BINARY_DEPTHTEST] = state;
	state = new OpenGLBinaryStateCullFace();		binaryState[IDX_BINARY_CULLFACE] = state;
	state = new OpenGLBinaryStateBlend();			binaryState[IDX_BINARY_BLEND] = state;
	state = new OpenGLBinaryStatePointSmooth();		binaryState[IDX_BINARY_POINTSMOOTH] = state;
	state = new OpenGLBinaryStateLineSmooth();		binaryState[IDX_BINARY_LINESMOOTH] = state;
	state = new OpenGLBinaryStatePolygonSmooth();	binaryState[IDX_BINARY_POLYGONSMOOTH] = state;
	state = new OpenGLBinaryStateColorMaterial();	binaryState[IDX_BINARY_COLORMATERIAL] = state;

	twoInputState.resize(numTwoInputState);

	OpenGLTwoInputStateVariable *state2 = NULL;

	state2 = new OpenGLTwoInputStateColorMaterial();	twoInputState[IDX_TWOINPUT_COLORMATERIAL] = state2;
	state2 = new OpenGLTwoInputStatePolygonMode();		twoInputState[IDX_TWOINPUT_POLYGONMODE] = state2;
	state2 = new OpenGLTwoInputStateBlendFunc();		twoInputState[IDX_TWOINPUT_BLENDFUNC] = state2;
}

void OpenGLState::glBind() const {
	for(int i=0;i<numBinaryState;++i)
		binaryState[i]->glBind();

	for(int i=0;i<numTwoInputState;++i)
		twoInputState[i]->glBind();
}