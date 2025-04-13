#include <stdafx.h>
#include <Nature/MassSpring/TestFunctions.h>
#include <Nature/MassSpring/Environment/Environments.h>
#include <Nature/MassSpring/Collision/Bridson02_SelfCollision.h>

using namespace std;

namespace springspace {

void test() {
	RandomGen<real> rnd;

	real t[3], w[2];

	for(size_t i=0;i<5;++i)
		Bridson02_SelfCollision::time_Coplanar(t, Vector3<real>(rnd), Vector3<real>(rnd), Vector3<real>(rnd), Vector3<real>(rnd), Vector3<real>(rnd), Vector3<real>(rnd));

	for(size_t i=0;i<5;++i)
		Bridson02_SelfCollision::project_vert2face(w, Vector3<real>(rnd), Vector3<real>(rnd), Vector3<real>(rnd), Vector3<real>(rnd));

	for(size_t i=0;i<20;++i)
		Bridson02_SelfCollision::closest_edge2edge(w, Vector3<real>(rnd), Vector3<real>(rnd), Vector3<real>(rnd), Vector3<real>(rnd));
	Bridson02_SelfCollision::closest_edge2edge(w, Vector3<real>(-1,0,0), Vector3<real>(1,0,0), Vector3<real>(0,-1,1), Vector3<real>(0,1,1));
	Bridson02_SelfCollision::closest_edge2edge(w, Vector3<real>(-1,0,0), Vector3<real>(1,0,0), Vector3<real>(-1,0,1), Vector3<real>(3,0,1));
}

};		// springspace