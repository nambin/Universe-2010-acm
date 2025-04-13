#include <stdafx.h>
#include <Nature/MassSpring/MassSpringIntegrator/MassSpringIntegrator.h>

using namespace springspace;

MassSpringIntegrator::MassSpringIntegrator() {
	initSimulation();
	set(1);
}
MassSpringIntegrator::MassSpringIntegrator(real _cfl) {
	initSimulation();
	set(_cfl);
}

void MassSpringIntegrator::initSimulation() {
	g = -9.81f;
}

void MassSpringIntegrator::set(real _cfl) {		
	maxCFL = _cfl;	
}

const real MassSpringIntegrator::getCFL(real dt, const QuadDynamicMesh<real>&) const {
	return maxCFL;
}

const real MassSpringIntegrator::getCFL(real dt, const TriDynamicMesh<real>&) const {
	return maxCFL;
}
