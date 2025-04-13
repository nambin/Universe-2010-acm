#include <stdafx.h>
#include <Principle/VirtualSpace/DrawLibrary/glDrawLibrary.h>

using namespace virtualspace;

void GlDrawLibrary::glDrawAxis(GLfloat endcoord, GLfloat width, const RGBA& cx, const RGBA& cy, const RGBA& cz) {

	glLineWidth(width);

	glColor4f(cx[0], cx[1], cx[2], cx[3]);
	glBegin(GL_LINES);
		glVertex3f(0, 0, 0);
		glVertex3f(endcoord, 0, 0);
	glEnd();

	glColor4f(cy[0], cy[1], cy[2], cy[3]);
	glBegin(GL_LINES);
		glVertex3f(0, 0, 0);
		glVertex3f(0, endcoord, 0);
	glEnd();

	glColor4f(cz[0], cz[1], cz[2], cz[3]);
	glBegin(GL_LINES);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 0, endcoord);
	glEnd();

}

void GlDrawLibrary::glDrawCube(GLfloat lx, GLfloat ly, GLfloat lz, const RGBA& cx, const RGBA& cy, const RGBA& cz) {

	glColor4f(cx[0], cx[1], cx[2], cx[3]);
	glBegin(GL_QUADS);
		glNormal3f(1, 0, 0);
		glVertex3f(lx,  0,  0);
		glVertex3f(lx, ly,  0);		
		glVertex3f(lx, ly, lz);
		glVertex3f(lx,  0, lz);		
		glNormal3f(-1, 0, 0);
		glVertex3f(0,  0,  0);
		glVertex3f(0,  0, lz);
		glVertex3f(0, ly, lz);
		glVertex3f(0, ly,  0);		
	glEnd();

	glColor4f(cy[0], cy[1], cy[2], cy[3]);
	glBegin(GL_QUADS);
		glNormal3f(0, 1, 0);
		glVertex3f( 0, ly,  0);
		glVertex3f( 0, ly, lz);
		glVertex3f(lx, ly, lz);
		glVertex3f(lx, ly,  0);		
		glNormal3f(0, -1, 0);
		glVertex3f(0,  0,  0);
		glVertex3f(lx, 0,  0);
		glVertex3f(lx, 0, lz);
		glVertex3f(0,  0, lz);						
	glEnd();

	glColor4f(cz[0], cz[1], cz[2], cz[3]);
	glBegin(GL_QUADS);
		glNormal3f(0, 0, 1);
		glVertex3f( 0,  0, lz);
		glVertex3f(lx,  0, lz);
		glVertex3f(lx, ly, lz);
		glVertex3f( 0, ly, lz);				
		glNormal3f(-1, 0, 0);
		glVertex3f( 0,  0,  0);
		glVertex3f( 0, ly,  0);
		glVertex3f(lx, ly,  0);
		glVertex3f(lx,  0,  0);			
	glEnd();

}

void GlDrawLibrary::glDrawCube_LineStrip(GLfloat Lx, GLfloat Ly, GLfloat Lz, const RGBA& color, GLfloat width) {
	glEnable(GL_COLOR_MATERIAL);

	glLineWidth(width);
	glColor4f( color[0], color[1], color[2], color[3] );

	glNormal3f(0.0f, 1.0f, 0.0f);

	glBegin( GL_LINES );
	glVertex3f( (GLfloat)0, (GLfloat)0, (GLfloat)0);
	glVertex3f( (GLfloat)Lx, (GLfloat)0, (GLfloat)0);
	glVertex3f( (GLfloat)0, (GLfloat)Ly, (GLfloat)0);
	glVertex3f( (GLfloat)Lx, (GLfloat)Ly, (GLfloat)0);
	glVertex3f( (GLfloat)0, (GLfloat)Ly, (GLfloat)Lz);
	glVertex3f( (GLfloat)Lx, (GLfloat)Ly, (GLfloat)Lz);
	glVertex3f( (GLfloat)0, (GLfloat)0, (GLfloat)Lz);
	glVertex3f( (GLfloat)Lx, (GLfloat)0, (GLfloat)Lz);
	glEnd();

	glBegin( GL_LINES );
	glVertex3f( (GLfloat)0, (GLfloat)0, (GLfloat)0);
	glVertex3f( (GLfloat)0, (GLfloat)Ly, (GLfloat)0);
	glVertex3f( (GLfloat)Lx, (GLfloat)0, (GLfloat)0);
	glVertex3f( (GLfloat)Lx, (GLfloat)Ly, (GLfloat)0);
	glVertex3f( (GLfloat)Lx, (GLfloat)0, (GLfloat)Lz);
	glVertex3f( (GLfloat)Lx, (GLfloat)Ly, (GLfloat)Lz);
	glVertex3f( (GLfloat)0, (GLfloat)0, (GLfloat)Lz);
	glVertex3f( (GLfloat)0, (GLfloat)Ly, (GLfloat)Lz);
	glEnd();

	glBegin( GL_LINES );
	glVertex3f( (GLfloat)0, (GLfloat)0, (GLfloat)0);
	glVertex3f( (GLfloat)0, (GLfloat)0, (GLfloat)Lz);
	glVertex3f( (GLfloat)Lx, (GLfloat)0, (GLfloat)0);
	glVertex3f( (GLfloat)Lx, (GLfloat)0, (GLfloat)Lz);
	glVertex3f( (GLfloat)0, (GLfloat)Ly, (GLfloat)0);
	glVertex3f( (GLfloat)0, (GLfloat)Ly, (GLfloat)Lz);
	glVertex3f( (GLfloat)Lx, (GLfloat)Ly, (GLfloat)0);
	glVertex3f( (GLfloat)Lx, (GLfloat)Ly, (GLfloat)Lz);		
	glEnd();

	glDisable(GL_COLOR_MATERIAL);
}

void GlDrawLibrary::glDrawViewOrientation(GLfloat windowWidth, GLfloat windowHeight, GLfloat ratio, const float* const mat) {

	glDisable(GL_LIGHTING);
	glDisable(GL_COLOR_MATERIAL);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_DEPTH_TEST);
	
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);

	const GLfloat _size = 1.0f;

	float _mat[16];
	memcpy(_mat, mat, sizeof(float) * 16);
	_mat[12] = 0.0f; _mat[13] = 0.0f; _mat[14] = -_size*3.0f;

	glViewport(0, 0, (GLsizei)windowWidth*ratio, (GLsizei)windowHeight*ratio);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glMultMatrixf(_mat);
	glTranslatef(-_size*0.5f, -_size*0.5f, -_size*0.5f);

	glDrawCube(_size, _size, _size);

	glPopMatrix();
}

void GlDrawLibrary::glDrawViewOrientation(GLfloat windowWidth, GLfloat windowHeight, GLfloat ratio, const double* const mat) {

	float _mat[16];
	for(int i=0;i<16;++i) _mat[i] = (float)mat[i];
	glDrawViewOrientation(windowWidth, windowHeight, ratio, _mat);

}