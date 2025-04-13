#include <stdafx.h>
#include "TestOpenGLScene.h"
#include <Principle/Geometry/Surface/OBJStruct.h>

using namespace std;
using namespace geometryspace;

void TestOpenGLScene::glBind() const {
	SceneObject::glBind();

	OpenGLState state = openGLState;
	state.glBind();

	glDraw();	
}

void TestOpenGLScene::glDraw() const {	

	float ratio = 1.0f;
	Light light1;	
	light1.setDiffuse(RGBA::white() * 1.3f * ratio);
	light1.setSpecular(RGBA::white() * 4.0f * ratio);
	light1.setPosition(Vector4f(0.0f, 0.5f, 0.0f, 1.0f));

	Light light2;	
	light2.setDiffuse(RGBA::white() * 1.3f * ratio);	
	light2.setPosition(Vector4f(0.0f, -0.5f, 0.0f, 1.0f));

	Light::initLights();
	light1.glBindLight();
	light2.glBindLight();

	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glEnable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glFrontFace(GL_CCW);

	glPolygonMode(GL_FRONT, GL_FILL);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

	glDisable(GL_COLOR_MATERIAL);

	//glFrontFace(GL_CW);

	//glShadeModel(GL_SMOOTH);
	//glEnable(GL_LIGHTING);
	//glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
	//glDisable(GL_CULL_FACE);
	//glDisable(GL_COLOR_MATERIAL);

//	Material material;
//	material.setDiffuse(RGBA::yellow() * 0.3f);
//	material.setSpecular(RGBA::yellow() * 0.3f);
//	material.glBind_FrontMaterial();
//
////	material.glBindMaterial();
//
//	material.setDiffuse(RGBA::blue() * 0.3f);
//	material.setSpecular(RGBA::blue() * 0.3f);
//	material.glBind_BackMaterial();
//
//	Vector3f v[3], n;
//	v[0] = Vector3f(1,0,0);
//	v[1] = Vector3f(0,0,-1);
//	v[2] = Vector3f(0,0,0);
//
////	glScalef(3.0f, 3.0f, 3.0f);
//
//	n = Vector3f(0,1,0);
//
//	glBegin(GL_TRIANGLES);
//	glNormal3fv(n.getptr());
//	for(int k=0;k<3;++k) 
//		glVertex3fv(v[k].getptr());
//	glEnd();

/*	n *= (-1.0f);
	glBegin(GL_TRIANGLES);
	glNormal3fv(n.getptr());
	for(int k=2;k>=0;--k) 
		glVertex3fv(v[k].getptr());
	glEnd();*/	
}

void TestOpenGLScene::resizeScreen(unsigned int width, unsigned int height) {
}

//////////////////////////////////////////////////////////////////////////
// Initial View Setting
const Vector3f TestOpenGLScene::viewEye() const { 
	return Vector3f(0,0,1);
}

const Vector3f TestOpenGLScene::viewRef() const { 
	return Vector3f(0,0,0);
}

const Vector3f TestOpenGLScene::viewUp() const  { 
	return Vector3f(0,1,0);
}

const float TestOpenGLScene::viewTheta() const  { 
	return 45;
}

TestOpenGLScene_IOHandler::TestOpenGLScene_IOHandler(TestOpenGLScene *_s) {
	sim = _s;
}

TestOpenGLScene_IOHandler::~TestOpenGLScene_IOHandler() {
}

void TestOpenGLScene_IOHandler::keyboardHandler(unsigned char key, int x, int y) {
	switch(key) {
		case 'o' : scene->getCamera()->set(sim->viewEye(), sim->viewRef(), sim->viewUp(), sim->viewTheta()); break;
		case '0' : scene->getCamera()->set(sim->viewEye_Side(), sim->viewRef_Side(), sim->viewUp_Side(), sim->viewTheta_Side()); break;

		case 'c' : scene->getCamera()->printCamera(); break;
	}
}

void TestOpenGLScene_IOHandler::mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y) {
}

void TestOpenGLScene_IOHandler::mouseMotionHandler(int x, int y) {
}
