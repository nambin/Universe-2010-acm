#include <stdafx.h>
#include "Header/TestFunction.h"
#include "Header/TestFunctionTmp.h"
#include "Header/Controller.h"

using namespace basis;
using namespace arithmetic;
using namespace virtualspace;
using namespace controlspace;

Controller *controller = NULL;
unsigned int millis = 1;

void initialize() {
	controller = new Controller();
	controller->set();
}

void render()	{
	controller->render();
}

void keyboardHandler(unsigned char key, int x, int y){
	controller->keyboardHandler(key, x, y);
}

void mouseHandler(int button, int state, int x, int y){
	controller->mouseHandler(button, state, x, y);
}
 
void mouseMotionHandler(int x, int y){
	controller->mouseMotionHandler(x,y);
}

void idleFunction() {
	controller->idleFunction();	
}

void timerFunction(int t) {
	controller->idleFunction();
	glutTimerFunc(millis, timerFunction, t);
}

void changeSize(GLsizei _width, GLsizei _height){
	controller->changeSize(_width, _height);
}

void main() {

	//testFunction();
	//arithmetic::CubicEquationSolverTest();	
	//springspace::test();

	// GLUT setting
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );	

	glutInitWindowSize(Scene::initWindowWidth(), Scene::initWindowHeight());
	glutInitWindowPosition(Scene::InitWindowX(), Scene::InitWindowY());
	glutCreateWindow("Universe");	

	glutDisplayFunc(render);	
	glutReshapeFunc(changeSize);		
	glutIdleFunc(idleFunction);
	//glutTimerFunc(millis, timerFunction, 1);

	glutKeyboardFunc(keyboardHandler);
	glutMouseFunc(mouseHandler);
	glutMotionFunc(mouseMotionHandler);						

	initialize();

	glutMainLoop();

}
