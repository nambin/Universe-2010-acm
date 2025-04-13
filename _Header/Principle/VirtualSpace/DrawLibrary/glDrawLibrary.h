#ifndef _VIRTUALSPACE_GL_DRAWLIBRARY_H_
#define _VIRTUALSPACE_GL_DRAWLIBRARY_H_

#include <stdlib.h>
#ifdef _WIN32
#include "gl/glut.h"
#else
#include <GLUT/glut.h>
#endif

#include "../../../Basis/Color/RGBA.h"
#include "../../Principle.h"

namespace virtualspace {

using basis::RGBA;

class PRINCIPLEDLL GlDrawLibrary {
public :
	static void glDrawAxis(GLfloat endcoord, GLfloat width=1.5f, const RGBA& cx=RGBA::red(), const RGBA& cy=RGBA::green(), const RGBA& cz=RGBA::blue());
	static void glDrawCube(GLfloat lx, GLfloat ly, GLfloat lz, const RGBA& cx=RGBA::red(), const RGBA& cy=RGBA::green(), const RGBA& cz=RGBA::blue());
	static void glDrawCube_LineStrip(GLfloat Lx, GLfloat Ly, GLfloat Lz, const RGBA& color = RGBA::white(), GLfloat width = 1.5f);

	static void glDrawViewOrientation(GLfloat windowWidth, GLfloat windowHeight, GLfloat ratio, const float* const mat);
	static void glDrawViewOrientation(GLfloat windowWidth, GLfloat windowHeight, GLfloat ratio, const double* const mat);

};
	
};	// namespace virtualspace

#endif