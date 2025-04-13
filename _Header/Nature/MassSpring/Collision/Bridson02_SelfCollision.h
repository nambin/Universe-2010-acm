#ifndef _MASS_SPRING_SPACE_BRIDSON02_SELF_COLLISION_H_
#define _MASS_SPRING_SPACE_BRIDSON02_SELF_COLLISION_H_

#include "../Collision/ClothSelfCollision.h"

namespace springspace {

class NATUREDLL Bridson02_SelfCollision : public ClothSelfCollision {
public : 
	struct Vert2Face {
		size_t	vert, face;
		real	w[2];			// barycentric coordinate
	};
	struct Edge2Edge {
		size_t	e0, e1;
		real	w[2];			
	};

public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Bridson02_SelfCollision();	
	~Bridson02_SelfCollision();

	//////////////////////////////////////////////////////////////////////////
	// Update	
	void treatCollision(real dt, QuadDynamicMesh<real>&);
	void treatCollision(real dt, TriDynamicMesh<real>&);	

	//////////////////////////////////////////////////////////////////////////
	// Resolve Proximity & Collision
	void resolve_Proximity(Vec3Array<real>& v, const Vec3Array<real>& x0, const vector<unsigned int>& v2f, const vector<unsigned int>& v2e, const Array<real>& m, const vector<Vector3<real>>& faceNormals) const;
	void resolve_GeometricCollision(Vec3Array<real>& v, const Vec3Array<real>& x0, const vector<unsigned int>& v2f, const vector<unsigned int>& v2e, const Array<real>& m, const vector<Vector3<real>>& faceNormals) const;

	//////////////////////////////////////////////////////////////////////////
	// Check Proximity & Collision
	void check_Proximity(vector<Vert2Face>&, vector<Edge2Edge>&, const BVH<unsigned int,real>*, const Vec3Array<real>& x, const vector<unsigned int>& v2f, const vector<unsigned int>& v2e) const;
	void check_Collision(vector<Vert2Face>&, vector<Edge2Edge>&, const BVH<unsigned int,real>*, const Vec3Array<real>& x, const Vec3Array<real>& v, const vector<unsigned int>& v2f, const vector<unsigned int>& v2e) const;

	static void project_vert2face(real w[2], const Vector3<real>& v, const Vector3<real>&, const Vector3<real>&, const Vector3<real>&); // w[2] is return value
	static void closest_edge2edge(real w[2], const Vector3<real>&  , const Vector3<real>&, const Vector3<real>&, const Vector3<real>&); // w[2] is return value
	static void closest_vert2edge(real& w  , const Vector3<real>& v, const Vector3<real>& e0, const Vector3<real>& dir);

	static void time_Coplanar(real t[3], const Vector3<real>&, const Vector3<real>&, const Vector3<real>&, const Vector3<real>&, const Vector3<real>&, const Vector3<real>&);

	//////////////////////////////////////////////////////////////////////////
	// Repulsion
	void get_repulsion(Vec3Array<real>& impulses, const vector<Vert2Face>&, const vector<Edge2Edge>&, const Vec3Array<real>& x, const Vec3Array<real>& v, const Array<real>& m, const vector<Vector3<real>>& faceNormals) const;

	real repulsion(const Vector3<real>& xn, const Vector3<real>& vn, const real m0, const real m1) const;
	real repulsion_inelastic(const Vector3<real>& vn, const real m0, const real m1) const;
	real repulsion_spring(const Vector3<real>& xn, const Vector3<real>& vn, const real m0, const real m1) const;
	real repulsion_friction(const Vector3<real>& vn_delta, const Vector3<real>& vt, const real mu, const real m0, const real m1) const;

	void interpolate_Impulse(Vec3Array<real>& impulses, const real impulse, const Vert2Face&, const vector<Vector3<real>>& faceNormals, const Array<real>& m) const;
	void interpolate_Impulse(Vec3Array<real>& impulses, const real impulse, const Edge2Edge&, const Array<real>& m) const;

	//////////////////////////////////////////////////////////////////////////
	// Rigid Impact Zone
	void apply_Rigid_Impact_Zone(Vec3Array<real>& v, const vector<Vert2Face>&, const vector<Edge2Edge>&, const Vec3Array<real>& x, const vector<unsigned int>& v2f, const vector<unsigned int>& v2e) const;

	//////////////////////////////////////////////////////////////////////////
	// Make BVH
	static BVH<unsigned int, real>* makeStaticBVH(const Vec3Array<real>& pos, const Array<pair<unsigned int, unsigned int>>&, const vector<unsigned int>& v2f, real eps = 1.0e-3);
	static BVH<unsigned int, real>* makeDynamicBVH(const Vec3Array<real>& pos, const Vec3Array<real>& vel, const Array<pair<unsigned int, unsigned int>>&, const vector<unsigned int>& v2f, real eps = 1.0e-6);
	static BVH<unsigned int, real>* mergeBVH(const Array<pair<unsigned int, unsigned int>>&, Array<BVH<unsigned int, real>*>&);

public : 
	real	thickness;		// cloth thickness
	real	k;				// spring-based repulsion

	mutable real dt;		// integration time
};
	
};		// springspace

#endif