#include <stdafx.h>
#include <Nature/MassSpring/Collision/Bridson02_SelfCollision.h>
#include <Principle/Arithmetic/ArrayUtil/ArrayUtil2.h>
#include <Principle/Arithmetic/Arithmetic/CubicEquationSolver_Cardano.h>

using namespace std;
using namespace arithmetic;
using namespace springspace;

// See [R.Bridson et al. 02]
void Bridson02_SelfCollision::check_Proximity(vector<Vert2Face>& vf, vector<Edge2Edge>& ee, const BVH<unsigned int,real>* bvh, const Vec3Array<real>& x, const vector<unsigned int>& v2f, const vector<unsigned int>& v2e) const {

}

void Bridson02_SelfCollision::check_Collision(vector<Vert2Face>& vf, vector<Edge2Edge>& ee, const BVH<unsigned int,real>* bvh, const Vec3Array<real>& x, const Vec3Array<real>& v, const vector<unsigned int>& v2f, const vector<unsigned int>& v2e) const {

}

void Bridson02_SelfCollision::project_vert2face(real w[2], const Vector3<real>& x4, const Vector3<real>& x1, const Vector3<real>& x2, const Vector3<real>& x3) {
	real eps = 1.0e-8;
	real A[2][2], b[2];

	Vector3<real> x13 = x1 - x3;
	Vector3<real> x23 = x2 - x3;
	Vector3<real> x43 = x4 - x3;

	A[0][0] = x13.sqrmag();
	A[1][1] = x23.sqrmag();
	A[0][1] = A[1][0] = dot(x13, x23);

	b[0] = dot(x13, x43);
	b[1] = dot(x23, x43);

	LinSol2x2(w, A, b, eps);

	// Verification _ It is verified
	//real w3 = 1 - w[0] - w[1];
	//Vector3<real> projected = w[0] * x1 + w[1] * x2 + w3 * x3;
	//Vector3<real> e0 = x2 - x1, e1 = x3 - x2, e2 = x1 - x3;
	//Vector3<real> l = projected - x4;
	//cout << dot(l, e0) << " " << dot(l, e1) << " " << dot(l, e2) << endl;
}

void Bridson02_SelfCollision::closest_edge2edge(real w[2], const Vector3<real>& x1, const Vector3<real>& x2, const Vector3<real>& x3, const Vector3<real>& x4) {
	real eps = 1.0e-8;

	real A[2][2], b[2];	

	Vector3<real> x21 = x2 - x1;	
	Vector3<real> x43 = x4 - x3;	
	Vector3<real> x31 = x3 - x1;

	if(sqrmag(cross(x21, x43)) < eps) {			
		real k = x43.mag() / x21.mag();
		real x = dot(x21, x31) / x21.sqrmag();
		real _min = MAX(real(0), x);
		real _max = MIN(real(1), x+k);
		w[0] = (_min + _max) * 0.5f;
		w[1] = (w[0] - x) / k;

		goto done;
	}

	A[0][0] = x21.sqrmag();
	A[1][1] = x43.sqrmag();
	A[0][1] = A[1][0] = - dot(x21, x43);

	b[0] = dot(x21, x31);
	b[1] = - dot(x43, x31);

	LinSol2x2(w, A, b, eps);

	if(0 <= w[0] && w[0] <= 1 && 0 <= w[1] && w[1] <= 1) 
		goto done;

	real moved[2];
	if(w[0] < 0)		{ moved[0] = w[0];		w[0] = 0; }
	else if(w[0] > 1)	{ moved[0] = w[0]-1;	w[0] = 1; }
	if(w[1] < 0)		{ moved[1] = w[1];		w[1] = 0; }
	else if(w[1] > 1)	{ moved[1] = w[1]-1;	w[1] = 1; }

	real movingLengthSqr[2];
	movingLengthSqr[0] = moved[0] * moved[0] * x21.sqrmag();
	movingLengthSqr[1] = moved[1] * moved[1] * x43.sqrmag();

	if(movingLengthSqr[0] < movingLengthSqr[1]) closest_vert2edge(w[1], x1 + w[0] * x21, x3, x4);
	else										closest_vert2edge(w[0], x3 + w[0] * x43, x1, x2);

done :
	// Verification _ It is verified
	//Vector3<real> v0 = x1 + w[0] * x21;
	//Vector3<real> v1 = x3 + w[1] * x43;
	//Vector3<real> l = v1 - v0;
	//cout << w[0] << " " << w[1] << " " << dot(l, x21) << " " << dot(l,x43) << endl;

	return;
}

void Bridson02_SelfCollision::closest_vert2edge(real& w, const Vector3<real>& v, const Vector3<real>& e0, const Vector3<real>& dir) {
	w = dot(v - e0, dir) / dir.sqrmag();
	if(w < 0)		w = 0;
	else if(w > 1)	w = 1;
}

void Bridson02_SelfCollision::time_Coplanar(real t[3], const Vector3<real>& x0, const Vector3<real>& v0, const Vector3<real>& x1, const Vector3<real>& v1, const Vector3<real>& x2, const Vector3<real>& v2) {

	static CubicEquationSolver_Cardano<real> cubicSolver;

	real a=0,b=0,c=0,d=0;
	real p=0,q=0,r=0;

	Vector3<real>::index_type j=0,k=0;
	for(Vector3<real>::index_type i=0;i<3;++i) {
		j = (i+1)%3; k = (i+2)%3;

		p = v0[j] * v1[k] - v1[j] * v0[k];
		q = x0[j] * v1[k] + x1[k] * v0[j] - x1[j] * v0[k] - x0[k] * v1[j];
		r = x0[j] * x1[k] - x0[k] * x1[j];

		a += p * v2[i];
		b += p * x2[i] + q * v2[i];
		c += q * x2[i] + r * v2[i];
		d += r * x2[i];
	}

	cubicSolver.solve(a,b,c,d);
	
	if(cubicSolver.numReal == 1)	  { t[0] = cubicSolver.sol[0].real();																		}
	else if(cubicSolver.numReal == 3) { t[0] = cubicSolver.sol[0].real(); t[1] = cubicSolver.sol[1].real(); t[2] = cubicSolver.sol[2].real();	}

	// Verification _ It is verified
	//cout << cubicSolver.numReal << endl;
	//for(size_t i=0;i<size_t(cubicSolver.numReal);++i) {
	//	real time = t[i];

	//	Vector3<real> k0 = x0 + time * v0;
	//	Vector3<real> k1 = x1 + time * v1;
	//	Vector3<real> k2 = x2 + time * v2;

	//	cout << time << " _ " << dot(cross(k0, k1), k2) << endl;
	//}
}
