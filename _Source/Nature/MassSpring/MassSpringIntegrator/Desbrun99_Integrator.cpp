#include <stdafx.h>
#include <Nature/MassSpring/MassSpringIntegrator/Desbrun99_Integrator.h>

using namespace std;
using namespace springspace;

Desbrun99_Integrator::Desbrun99_Integrator() { init(); }
Desbrun99_Integrator::Desbrun99_Integrator(real _cfl) 
: MassSpringIntegrator(_cfl) { init(); }

Desbrun99_Integrator::~Desbrun99_Integrator() {}

void Desbrun99_Integrator::init() {
	MassSpringIntegrator::initSimulation();
	filter.init();	
}

const real Desbrun99_Integrator::getCFL(real dt, const QuadDynamicMesh<real>&) const { return maxCFL; }
const real Desbrun99_Integrator::getCFL(real dt, const TriDynamicMesh<real>&) const	 { return maxCFL; }

void Desbrun99_Integrator::update(real dt, QuadDynamicMesh<real>& mesh) {
	if(filter.dimRow() == 0 || filter.dimCol() == 0)
		basis::throwError("Error in Desbrun99_Integrator::update() _ Should make filter first\n");
	
	Vec3Array<real> force(mesh.size());
	Vector3<real> com = mesh.centerOfMass();

	swap(mesh.pos, mesh.pos0);
	swap(mesh.vel, mesh.vel0);

	addSpringForce(force, dt, mesh);
	filterForce(force);

	Vector3<real> torque = getAdditionalGlobalTorque(force, mesh.pos0);	

	cout << "Torque : " << torque << endl;

	updateDynamicStates(mesh, dt, force);

	printMaxVector3(mesh.pos);

//	postCorrectionAngularMomentum(mesh, dt, torque, com);

	printMaxVector3(mesh.pos);	

//	postStep_IK(mesh, dt);	//updateFinalVelocity(mesh, dt);	
}

void Desbrun99_Integrator::update(real dt, TriDynamicMesh<real>& mesh) {
	basis::throwError("Not Implemented _ Desbrun99_Integrator::update()\n");
}
