#ifndef _VIRTUALSPACE_OPENGL_STATE_VARIABLE_H_
#define _VIRTUALSPACE_OPENGL_STATE_VARIABLE_H_

#include <stdlib.h>
#ifdef _WIN32
#include "gl/glut.h"
#else
#include <GLUT/glut.h>
#endif

namespace virtualspace {

class OpenGLStateVariable {
	virtual void glBind() const = 0;
};

class OpenGLTwoInputStateVariable : public OpenGLStateVariable {
public :
	void set(GLenum _v0, GLenum _v1) {
		var0 = _v0; var1 = _v1;
	}
	virtual void glBind() const = 0;

	GLenum var0, var1;
};

class OpenGLTwoInputStateColorMaterial : public OpenGLTwoInputStateVariable {
public : 
	OpenGLTwoInputStateColorMaterial() { var0 = GL_FRONT_AND_BACK; var1 = GL_AMBIENT_AND_DIFFUSE; }

	void glBind() const {
		glColorMaterial(var0, var1);
	}
};

class OpenGLTwoInputStatePolygonMode : public OpenGLTwoInputStateVariable {
public : 
	OpenGLTwoInputStatePolygonMode() { var0 = GL_FRONT_AND_BACK; var1 = GL_FILL; }

	void glBind() const {
		glPolygonMode(var0, var1);
	}
};

class OpenGLTwoInputStateBlendFunc : public OpenGLTwoInputStateVariable {
public : 
	OpenGLTwoInputStateBlendFunc() { var0 = GL_SRC_ALPHA; var1 = GL_ONE_MINUS_SRC_ALPHA; }

	void glBind() const {
		glBlendFunc(var0, var1);
	}
};

};			// namespace virtualspace;

#endif
