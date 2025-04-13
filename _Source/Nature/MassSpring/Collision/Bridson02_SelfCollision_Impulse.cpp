#include <stdafx.h>
#include <Nature/MassSpring/Collision/Bridson02_SelfCollision.h>

using namespace std;
using namespace springspace;

void Bridson02_SelfCollision::get_repulsion(Vec3Array<real>& impulses, const vector<Vert2Face>&, const vector<Edge2Edge>&, const Vec3Array<real>& x, const Vec3Array<real>& v, const Array<real>& m, const vector<Vector3<real>>& faceNormals) const {

}

real Bridson02_SelfCollision::repulsion(const Vector3<real>& xn, const Vector3<real>& vn, const real m0, const real m1) const {
	return 0;
}

real Bridson02_SelfCollision::repulsion_inelastic(const Vector3<real>& vn, const real m0, const real m1) const {
	return 0;
}

real Bridson02_SelfCollision::repulsion_spring(const Vector3<real>& xn, const Vector3<real>& vn,const real m0, const real m1) const {
	return 0;
}

real Bridson02_SelfCollision::repulsion_friction(const Vector3<real>& vn_delta, const Vector3<real>& vt, const real mu, const real m0, const real m1) const {
	return 0;
}

void Bridson02_SelfCollision::interpolate_Impulse(Vec3Array<real>& impulses, const real impulse, const Vert2Face&, const vector<Vector3<real>>& faceNormals, const Array<real>& mass) const {

}

void Bridson02_SelfCollision::interpolate_Impulse(Vec3Array<real>& impulses, const real impulse, const Edge2Edge&, const Array<real>& mass) const {

}