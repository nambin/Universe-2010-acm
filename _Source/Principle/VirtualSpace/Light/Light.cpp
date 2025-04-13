#include <stdafx.h>
#include <Principle/VIrtualSpace/Light/Light.h>

using namespace basis;
using namespace virtualspace;

void Light::glDraw(float size) const {
	
	glDisable(GL_COLOR_MATERIAL);
	glEnable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);

	GLUquadricObj *qobj = gluNewQuadric();

	gluQuadricDrawStyle(qobj, (GLenum) GLU_FILL);
	gluQuadricNormals(qobj, (GLenum) GLU_SMOOTH);

	GLfloat black[] = {0,0,0,1};
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, black);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, black);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, black);
	glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, RGBA::white().getptr());

	glTranslatef(position[0], position[1], position[2]);		
	gluSphere(qobj, size, 15, 15);
	glTranslatef(-position[0], -position[1], -position[2]);
}

int Light::lightCount = 0;
void Light::glBindLight() const {

	if(!enabled)
		return;

	glEnable(GL_LIGHTING);

	lightID = -1;
	switch(Light::lightCount)	{
		case 0: lightID = GL_LIGHT0; break;
		case 1: lightID = GL_LIGHT1; break;
		case 2: lightID = GL_LIGHT2; break;
		case 3: lightID = GL_LIGHT3; break;
		case 4: lightID = GL_LIGHT4; break;
		case 5: lightID = GL_LIGHT5; break;
		case 6: lightID = GL_LIGHT6; break;
		case 7: lightID = GL_LIGHT7; break;
		default:
			throwError("Can't use more than 8 lights");
	}

	if(lightID != -1)
		++Light::lightCount;
	
	glLightf(lightID, GL_CONSTANT_ATTENUATION , attenuation[0]);
	glLightf(lightID, GL_LINEAR_ATTENUATION   , attenuation[1]);
	glLightf(lightID, GL_QUADRATIC_ATTENUATION, attenuation[2]);

	glLightfv(lightID, GL_POSITION, position.getptr());
	glLightfv(lightID, GL_AMBIENT,  ambient.getptr());
	glLightfv(lightID, GL_DIFFUSE,  diffuse.getptr());
	glLightfv(lightID, GL_SPECULAR, specular.getptr());

	glEnable(lightID);
}

void Light::glUnbindLight() const {
	if(!enabled)
		return;

	glDisable(lightID);
}