#include <stdafx.h>
#include <Illumination/RealtimeShader/cgShader/cgShader_VolumeRendering.h>

using namespace std;
using namespace virtualspace;
using namespace shaderspace;

void cgShader_VolumeRendering::glDraw() const {
	glDisable(GL_LIGHTING);

	enableRenderbuffers();
	renderBackFace();
	disableRenderBuffers();

	rayCasting();	

	if(dBoundingBox)
		GlDrawLibrary::glDrawCube_LineStrip(Lx, Ly, Lz, RGBA::white());
}

void cgShader_VolumeRendering::enableRenderbuffers() const {
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, frameBuffer);
}

void cgShader_VolumeRendering::disableRenderBuffers() const {
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}

void cgShader_VolumeRendering::renderBackFace() const {
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, backFaceBuffer, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_FRONT);
	drawQuads(Lx, Ly, Lz);
	glDisable(GL_CULL_FACE);
}

void cgShader_VolumeRendering::rayCasting() const {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	cgGLEnableProfile(vertexProfile);
	cgGLEnableProfile(fragmentProfile);
	cgGLBindProgram(vertexProgram);
	cgGLBindProgram(fragmentProgram);
	cgGLSetParameter1f(cgGetNamedParameter( fragmentProgram, "stepsize") , stepsize);
	
	cgShaderLibrary::set_TextureParameter(param1, "tex", backFaceBuffer, fragmentProgram);
	cgShaderLibrary::set_TextureParameter(param2, "volume_tex", volume_texture, fragmentProgram);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	drawQuads(Lx, Ly, Lz);
	glDisable(GL_CULL_FACE);

	cgGLDisableProfile(vertexProfile);
	cgGLDisableProfile(fragmentProfile);
}

void cgShader_VolumeRendering::drawVertex(float x, float y, float z) const {
	float xc = x / Lx;
	float yc = y / Ly;
	float zc = z / Lz;

	glColor3f(xc, yc, zc);
	glMultiTexCoord3f(GL_TEXTURE1, xc, yc, zc);
	glVertex3f(x,y,z);
}

void cgShader_VolumeRendering::drawQuads(float x, float y, float z) const {

	glBegin(GL_QUADS);
		// Back side 
		glNormal3f(0, 0, -1);
		drawVertex(0, 0, 0);
		drawVertex(0, y, 0);
		drawVertex(x, y, 0);
		drawVertex(x, 0, 0);

		// Front side 
		glNormal3f(0, 0, 1);
		drawVertex(0, 0, z);
		drawVertex(x, 0, z);
		drawVertex(x, y, z);
		drawVertex(0, y, z);

		// Top side 
		glNormal3f(0, 1, 0);
		drawVertex(0, y, 0);
		drawVertex(0, y, z);
		drawVertex(x, y, z);
		drawVertex(x, y, 0);

		// Bottom side 
		glNormal3f(0, -1, 0);
		drawVertex(0, 0, 0);
		drawVertex(x, 0, 0);
		drawVertex(x, 0, z);
		drawVertex(0, 0, z);

		// Left side 
		glNormal3f(-1, 0, 0);
		drawVertex(0, 0, 0);
		drawVertex(0, 0, z);
		drawVertex(0, y, z);
		drawVertex(0, y, 0);

		// Right side 
		glNormal3f(1, 0, 0);
		drawVertex(x, 0, 0);
		drawVertex(x, y, 0);
		drawVertex(x, y, z);
		drawVertex(x, 0, z);
	glEnd();

}