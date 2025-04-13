#include <stdafx.h>
#include <Nature/MassSpring/ClothScene/ClothSceneIntegrator.h>

using namespace basis;
using namespace springspace;

ClothSceneIntegrator_Basic::ClothSceneIntegrator_Basic(real _dt) : ClothSceneIntegrator(_dt) {
}

void ClothSceneIntegrator_Basic::update(real dt, ClothScene_Objects* objects, DynamicMesh* mesh) {

	real dtSum = 0;
	bool isEnd = false;

	while(!isEnd) {
		real _dt = dt - dtSum;

		if(dt0 < _dt)	_dt = dt0;
		else			isEnd = true;
		
		if(dynamic_cast<TriDynamicMesh<real>*>(mesh) != NULL)	update(_dt, objects, reinterpret_cast<TriDynamicMesh<real>*>(mesh));
		if(dynamic_cast<QuadDynamicMesh<real>*>(mesh) != NULL)	update(_dt, objects, reinterpret_cast<QuadDynamicMesh<real>*>(mesh));

		dtSum += _dt;
	}
}

void ClothSceneIntegrator_Basic::update(real dt, ClothScene_Objects* objects, QuadDynamicMesh<real>* mesh) {
	Timer timer;

	real dtSum = 0;
	size_t nSubSteps = 0;

	while(true) {
		real _dt = dt - dtSum;
		real cfl = dynamics->getCFL(_dt, *mesh);
		nSubSteps++;

		if(cfl <= dynamics->maxCFL) {
			dynamics->update(_dt, *mesh);
			break;
		}
		else {
			_dt = _dt / cfl * dynamics->maxCFL;
			dynamics->update(_dt, *mesh);
			dtSum += _dt;
		}
	}	

	cout << "Strain : " << mesh->maxStrain() << endl;

	if(postDynamics) {
		postDynamics->postProcess(dt, *mesh);
		cout << "Strain : " << mesh->maxStrain() << endl;
	}

	if(objects && clothObjectCollision) {
		mesh->constraint.releaseSoftConstraints();
		clothObjectCollision->treatCollision(*objects, dt, *mesh);
	}	

	timer.printTime("QuadDynamicMeshSimulator::update");
}

void ClothSceneIntegrator_Basic::update(real dt, ClothScene_Objects* objects, TriDynamicMesh<real>* mesh) {
	Timer timer;

	real dtSum = 0;
	size_t nSubSteps = 0;

	while(true) {
		real _dt = dt - dtSum;
		real cfl = dynamics->getCFL(_dt, *mesh);
		nSubSteps++;

		if(cfl <= dynamics->maxCFL) {
			dynamics->update(_dt, *mesh);
			break;
		}
		else {
			_dt = _dt / cfl * dynamics->maxCFL;
			dynamics->update(_dt, *mesh);
			dtSum += _dt;
		}
	}	

	if(postDynamics) {
		postDynamics->postProcess(dt, *mesh);
	}

	if(objects && clothObjectCollision) {
		mesh->constraint.releaseSoftConstraints();
		clothObjectCollision->treatCollision(*objects, dt, *mesh);
	}

	cout << "SubStep : " << nSubSteps;
	timer.printTime(" _ Spent");	
}