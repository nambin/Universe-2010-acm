#ifndef _VIRTUALSPACE_OPENGL_BINARYSTATE_VARIABLE_H_
#define _VIRTUALSPACE_OPENGL_BINARYSTATE_VARIABLE_H_

#include <stdlib.h>
#ifdef _WIN32
#include "gl/glut.h"
#else
#include <GLUT/glut.h>
#endif

#include "OpenGLStateVariable.h"

namespace virtualspace {

class OpenGLBinaryStateVariable : public OpenGLStateVariable {
public : 
	virtual void glBind() const = 0;

	bool enable;
};

class OpenGLBinaryStateLighting : public OpenGLBinaryStateVariable {
public : 
	OpenGLBinaryStateLighting() { enable = true; }

	void glBind() const {
		if(enable)	glEnable(GL_LIGHTING);
		else		glDisable(GL_LIGHTING);
	}
};

class OpenGLBinaryStateDepthTest : public OpenGLBinaryStateVariable {
public : 
	OpenGLBinaryStateDepthTest() { enable = true; }

	void glBind() const {
		if(enable)	glEnable(GL_DEPTH_TEST);
		else		glDisable(GL_DEPTH_TEST);
	}
};

class OpenGLBinaryStateCullFace : public OpenGLBinaryStateVariable {
public : 
	OpenGLBinaryStateCullFace() { enable = true; }

	void glBind() const {
		if(enable)	glEnable(GL_CULL_FACE);
		else		glDisable(GL_CULL_FACE);
	}
};

class OpenGLBinaryStateBlend : public OpenGLBinaryStateVariable {
public : 
	OpenGLBinaryStateBlend() { enable = true; }

	void glBind() const {
		if(enable)	glEnable(GL_BLEND);
		else		glDisable(GL_BLEND);
	}
};

class OpenGLBinaryStatePointSmooth : public OpenGLBinaryStateVariable {
public : 
	OpenGLBinaryStatePointSmooth() { enable = true; }

	void glBind() const {
		if(enable)	glEnable(GL_POINT_SMOOTH);
		else		glDisable(GL_POINT_SMOOTH);
	}
};

class OpenGLBinaryStateLineSmooth : public OpenGLBinaryStateVariable {
public : 
	OpenGLBinaryStateLineSmooth() { enable = true; }

	void glBind() const {
		if(enable)	glEnable(GL_LINE_SMOOTH);
		else		glDisable(GL_LINE_SMOOTH);
	}
};

class OpenGLBinaryStatePolygonSmooth : public OpenGLBinaryStateVariable {
public : 
	OpenGLBinaryStatePolygonSmooth() { enable = true; }

	void glBind() const {
		if(enable)	glEnable(GL_POLYGON_SMOOTH);
		else		glDisable(GL_POLYGON_SMOOTH);
	}
};

class OpenGLBinaryStateColorMaterial : public OpenGLBinaryStateVariable {
public : 
	OpenGLBinaryStateColorMaterial() { enable = false; }

	void glBind() const {
		if(enable)	glEnable(GL_COLOR_MATERIAL);
		else		glDisable(GL_COLOR_MATERIAL);
	}
};

};			// namespace virtualspace;

#endif
