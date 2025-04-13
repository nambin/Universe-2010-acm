#include <stdafx.h>
#include "Header/Controller.h"
#include "Header/SceneNodeAdder.h"

Controller::Controller()	{ init();		}
Controller::~Controller()	{ release();	}
void Controller::init()		{ scene = NULL;	}

void Controller::set() {

	scene = new Scene();

	InitLights(scene);
	SRL_Experiment(scene);
	
//	AddFluidRelated(scene);
//	AddMassSpringRelated(scene);
//	AddGeometryRelated(scene);

//	illuminationspace::test();

//	AddOpenGLTest(scene);

	cout << "controller" << endl;
}

void Controller::release() {
	if(scene)
		delete scene; 
	scene = NULL;
}

