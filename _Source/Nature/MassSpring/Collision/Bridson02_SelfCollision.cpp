#include <stdafx.h>
#include <Nature/MassSpring/Collision/Bridson02_SelfCollision.h>

using namespace std;
using namespace springspace;

Bridson02_SelfCollision::Bridson02_SelfCollision() {}
Bridson02_SelfCollision::~Bridson02_SelfCollision() {}

//////////////////////////////////////////////////////////////////////////
// Update	
void Bridson02_SelfCollision::treatCollision(real dt, QuadDynamicMesh<real>& mesh) {
	basis::throwError("Not Implemented _ Bridson02_SelfCollision::treatCollision()");
}

void Bridson02_SelfCollision::treatCollision(real dt, TriDynamicMesh<real>& mesh) {
}

void Bridson02_SelfCollision::resolve_Proximity(Vec3Array<real>& v, const Vec3Array<real>& x0, const vector<unsigned int>& v2f, const vector<unsigned int>& v2e, const Array<real>& m, const vector<Vector3<real>>& faceNormals) const {

}

void Bridson02_SelfCollision::resolve_GeometricCollision(Vec3Array<real>& v, const Vec3Array<real>& x0, const vector<unsigned int>& v2f, const vector<unsigned int>& v2e, const Array<real>& m, const vector<Vector3<real>>& faceNormals) const {

}