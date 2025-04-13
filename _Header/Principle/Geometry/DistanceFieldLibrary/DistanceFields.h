#ifndef _GEOMETRYSPACE_DISTANCE_FIELDS_H_
#define _GEOMETRYSPACE_DISTANCE_FIELDS_H_

#include "../../Principle.h"
#include "../../Arithmetic/Vector/Vector2.h"
#include "../../Arithmetic/Vector/Vector3.h"
#include "../../Geometry/Surface/TriMesh_Simple.h"
#include "../../Arithmetic/GridLibrary/Grid3.h"
//#include "../../Geometry/Levelset/Levelset_Redistancing.h"

namespace geometryspace {

using arithmetic::Vector2;
using arithmetic::Vector3;

template<typename T> 
class DistanceField {
public : 
	virtual ~DistanceField() {}

	const T operator()(const Vector2<T> &pos) const { return getPhi(pos); }
	const T operator()(const Vector3<T> &pos) const { return getPhi(pos); }

	virtual const T getPhi(const Vector2<T> &pos) const { throwError("Not implemented _ DistanceField::getPhi(const Vector2<T> &pos)"); return 1.0e+6; }
	virtual const T getPhi(const Vector3<T> &pos) const { throwError("Not implemented _ DistanceField::getPhi(const Vector3<T> &pos)"); return 1.0e+6; }

	virtual const unsigned int getDimension() const = 0;
};

template<typename T> 
class DistanceField2D : public DistanceField<T> {
public : 
	virtual ~DistanceField2D() {}	
	const unsigned int getDimension() const { return 2; }
};

template<typename T> 
class CircleDistanceField2D : public DistanceField2D<T> {
public : 
	CircleDistanceField2D(const Vector2<T> &_c, T _r) { init(_c, _r); }
	void init(const Vector2<T> &_c, T _r) { center = _c; radius = _r; }
	const T getPhi(const Vector2<T> &pos) const {
		return length(pos, center) - radius;
	}

private : 
	Vector2<T>	center;
	T			radius;
};

template<typename T> 
class RectangleDistanceField2D : public DistanceField2D<T> {
public : 
	RectangleDistanceField2D(const Vector2<T> &__min, const Vector2<T> &__max) { init(__min, __max); }
	void init(const Vector2<T> &__min, const Vector2<T> &__max) {
		_min = minVec(__min, __max);
		_max = maxVec(__min, __max);
	}
	const T getPhi(const Vector2<T> &pos) const {
		T sqDist = 0.0f;
		for(Vector2<T>::index_type i=0;i<2;++i) {
			T v = pos[i];
			if(v < _min[i]) sqDist += (_min[i] - v) * (_min[i] - v);
			if(v > _max[i]) sqDist += (v - _max[i]) * (v - _max[i]);
		}
		if(sqDist != 0.0f) return sqrt(sqDist);

		T _minDist = 1.0e+8;
		for(Vector2<T>::index_type i=0;i<2;++i) {
			T v = pos[i];
			_minDist = MIN(_minDist, v - _min[i]);
			_minDist = MIN(_minDist, _max[i] - v);
		}
		return -_minDist;
	}

private : 
	Vector2<T> _min, _max;
};

template<typename T> 
class ZalesakDistanceField2D : public DistanceField2D<T> {
public : 
	ZalesakDistanceField2D(T _L) : L(_L), circle(Vector2<T>(0), 1), rec(Vector2<T>(0), Vector2<T>(0)) { init(_L); }
	void init(T _L) {
		L = _L;
		circle.init(Vector2<T>(0.5f, 0.75f) * L, 0.15f*L);

		Vector2<T> _min, _max;
		_min[0] = 0.5f*L - 0.025f*L;
		_min[1] = 0.75f*L - sqrt(SQ(0.15f) - SQ(0.025f)) * L;
		_max = _min + Vector2<T>(0.05f, 0.25f) * L;

		rec.init(_min - Vector2<T>(0,0.1f*L), _max);
	}
	const T getPhi(const Vector2<T> &pos) const {
		T p0 = circle.getPhi(pos);
		T p1 = -rec.getPhi(pos);
		return MAX(p0, p1);
	}

private : 
	T								L;
	CircleDistanceField2D<T>		circle;
	RectangleDistanceField2D<T>		rec;
};

//template<typename T> 
//class ZalesakDistanceField2D_Arbitrary : public ZalesakDistanceField2D<T> {
//public : 
//	ZalesakDistanceField2D_Arbitrary(const Vector2<T>& c, const T r) : ZalesakDistanceField2D(1) {} 
//	
//	const T getPhi(const Vector2<T> &pos) const {
//		
//	}
//
//private : 
//	const Vector2<T>	center;
//	const T				r;	
//};

template<typename T> 
class RussoDistanceField2D : public DistanceField2D<T> {
public : 
	RussoDistanceField2D(const Vector2<T>& _c = Vector2<T>(0), const T s = T(1.0)) : center(_c), scale(s) {}	
	const T getPhi(const Vector2<T> &pos) const {
		T x = pos[0] - center[0], y = pos[1] - center[1];
		return ((0.1 + SQ(x-1) + SQ(y-1))) * (sqrt(x*x + y*y) - 1) * scale;
	}

private : 
	Vector2<T>	center;
	T			scale;
};

template<typename T> 
class WaterDropDistanceField2D : public DistanceField2D<T> {
public : 
	WaterDropDistanceField2D(const Vector2<T> &_c, T _r, T _h) : circle(_c, _r), height(_h) {}
	void init(const Vector2<T> &_c, T _r, T _h) {
		circle.init(_c, _r);
		height = _h;
	}
	const T getPhi(const Vector2<T> &pos) const {
		T phi0 = pos[1] - height;
		T phi1 = circle.getPhi(pos);

		return (abs(phi0) < abs(phi1)) ? phi0 : phi1;
	}

private : 
	CircleDistanceField2D<T>	circle;
	T							height;
};

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

template<typename T> 
class DistanceField3D : public DistanceField<T> {
public : 
	virtual ~DistanceField3D() {}
	const unsigned int getDimension() const { return 3; }
	
	//virtual const AABB<T> getAABB() const {
	//	return AABB<T>(Vector3<T>(-1.0e+8), Vector3<T>(1.0e+8));
	//}
};

template<typename T>
class SphereDistanceField3D : public DistanceField3D<T> {
public : 
	SphereDistanceField3D(const Vector3<T> &_c, T _r) { init(_c, _r); }
	void init(const Vector3<T> &_c, T _r) { center = _c; radius = _r; }
	const T getPhi(const Vector3<T> &pos) const {
		return length(pos, center) - radius;
	}

private : 
	Vector3<T>	center;
	T			radius;
};

template<typename T> 
class CubeDistanceField3D : public DistanceField3D<T> {
public : 
	CubeDistanceField3D(const Vector3<T> &__min, const Vector3<T> &__max) { init(__min, __max); }

	const Vector3<T> get_min() const { return _min; }
	const Vector3<T> get_max() const { return _max; }

	void init(const Vector3<T> &__min, const Vector3<T> &__max) {
		_min = minVec(__min, __max);
		_max = maxVec(__min, __max);
	}
	const T getPhi(const Vector3<T> &pos) const {
		T sqDist = 0.0f;
		for(unsigned int i=0;i<3;++i) {
			T v = pos[i];
			if(v < _min[i]) sqDist += (_min[i] - v) * (_min[i] - v);
			if(v > _max[i]) sqDist += (v - _max[i]) * (v - _max[i]);
		}
		if(sqDist != 0.0f) return sqrt(sqDist);

		T _minDist = 1.0e+8;
		for(unsigned int i=0;i<3;++i) {
			T v = pos[i];
			_minDist = MIN(_minDist, v - _min[i]);
			_minDist = MIN(_minDist, _max[i] - v);
		}
		return -_minDist;
	}

private : 
	Vector3<T> _min, _max;
};

template<typename T>
class CylinderDistanceField : public DistanceField3D<T> {		// center is (0,0,0), y-direction
public : 
	CylinderDistanceField(T r, T len) : radius(r), length_half(len*0.5) {}
	
	void set_rigid_transformation(const arithmetic::Transform<T>& transform) { rigid_transform = transform; }

	const T							get_radius() const { return radius;				}
	const T							get_length() const { return length_half * 2.0;	}
	const arithmetic::Transform<T>& get_rigid_transform() const { return rigid_transform; }

	const T getPhi(const Vector3<T>& pos) const {
		const Vector3<T> p = rigid_transform.getInverse().transformPoint(pos);
		if(-length_half < p[1] && p[1] < length_half)
			return sqrt(p[0]*p[0] + p[2]*p[2]) - radius;

		const T sqrmag = (p[0]*p[0] + p[2]*p[2]);
		if(sqrmag < radius * radius)
			return ::abs(p[1] - length_half * SIGN(p[1]));

		Vector3<T> closest; closest[1] = length_half * SIGN(p[1]);
		closest[0] = p[0] / sqrt(sqrmag);
		closest[2] = p[2] / sqrt(sqrmag);

		return length(closest, p);
	}

protected : 
	arithmetic::Transform<T>	rigid_transform;
	T							radius, length_half;
};

template<typename T> 
class ZalesakDistanceField3D : public DistanceField3D<T> {
public : 
	ZalesakDistanceField3D(T _L) : L(_L), sphere(Vector3<T>(0), 1), cube(Vector3<T>(0), Vector3<T>(0)) { init(_L); }
	void init(T _L) {
		L = _L;
		sphere.init(Vector3<T>(0.5, 0.75, 0.5) * L, 0.15 * L);

		Vector3<T> _min, _max;
		_min[0] = T(0.5 - 0.025) * L;
		_min[1] = T(0.75) * L - sqrt(SQ(0.15) - SQ(0.025)) * L;
		_min[2] = T(0.0);
		_max = _min + Vector3<T>(0.05, 0.25, 1.0) * L;

		cube.init(_min - Vector3<T>(0,0.1*L,0), _max);
	}
	const T getPhi(const Vector3<T> &pos) const {
		return MAX(sphere.getPhi(pos), - cube.getPhi(pos));
	}

private : 
	T							L;
	SphereDistanceField3D<T>	sphere;
	CubeDistanceField3D<T>		cube;
};

template<typename T> 
class ZalesakDistanceField3D_Enright : public DistanceField3D<T> {
public : 
	ZalesakDistanceField3D_Enright(T _L) : L(_L), sphere(Vector3<T>(0), 1), cube(Vector3<T>(0), Vector3<T>(0)) { init(_L); }
	void init(T _L) {
		L = _L;
		sphere.init(Vector3<T>(0.5, 0.75, 0.5) * L, 0.15 * L);

		Vector3<T> _min, _max;
		_min[0] = T(0.5 - 0.025) * L;
		_min[1] = T(0.75) * L - sqrt(SQ(0.15) - SQ(0.025)) * L;
		_min[2] = T(0.0);
		_max = _min + Vector3<T>(0.05, 0.125, 1.0) * L;

		cube.init(_min - Vector3<T>(0,0.1*L,0), _max);
	}
	const T getPhi(const Vector3<T> &pos) const {
		return MAX(sphere.getPhi(pos), - cube.getPhi(pos));
	}

private : 
	T							L;
	SphereDistanceField3D<T>	sphere;
	CubeDistanceField3D<T>		cube;
};

template<typename T> 
class WaterDropDistanceField3D : public DistanceField3D<T> {
public : 
	WaterDropDistanceField3D(const Vector3<T> &_c, T _r, T _h) : sphere(_c, _r), height(_h) {}
	void init(const Vector3<T> &_c, T _r, T _h) {
		sphere.init(_c, _r);
		height = _h;
	}
	const T getPhi(const Vector3<T> &pos) const {
		T phi0 = pos[1] - height;
		T phi1 = sphere.getPhi(pos);

		return (abs(phi0) < abs(phi1)) ? phi0 : phi1;
	}

private : 
	SphereDistanceField3D<T>	sphere;
	T							height;
};

template<typename T> 
class ObjectDistanceField3D : public DistanceField3D<T> {
public : 
	typedef TriMesh_Simple<Vector3<T>, vector<Vector3<T>>>	Mesh;

	ObjectDistanceField3D() : octree(NULL), scale_factor(1.0) {}
	ObjectDistanceField3D(const std::string& file, bool unitize = true, bool change_order = false) : octree(NULL), scale_factor(1.0) { init(file, unitize, change_order); }
	~ObjectDistanceField3D() { if(octree) delete octree; }

	const AABB<T> getAABB() const {
		AABB<T> aabb = mesh.getAABB(transform);
		aabb.expand(aabb.extent() * T(0.1));
		
		return aabb;
	}

	void set_scale(const T s)											 { scale_factor = s;		transform = rigid_transform * arithmetic::scale(scale_factor, scale_factor, scale_factor); }
	void set_rigid_transformation(const arithmetic::Transform<T>& trans) { rigid_transform = trans;	transform = rigid_transform * arithmetic::scale(scale_factor, scale_factor, scale_factor); }

	void init(const std::string& file, bool unitize, bool change_order) { 
		Timer timer;

		mesh.load(file.c_str());

		if(unitize)			mesh.unitize();
		if(change_order)	mesh.changeTriangleOrder();

		mesh.makeFaceNormals();
		mesh.makeVertexNormals();

		timer.logTime("ObjectDistanceField3D::init() _ mesh is loaded", log_info); 

		t2v.resize(mesh.numVertices);
		for(unsigned int i=0;i<mesh.numFaces;i++)	{
			t2v[mesh.v2f[i*3  ]].push_back(i);
			t2v[mesh.v2f[i*3+1]].push_back(i);
			t2v[mesh.v2f[i*3+2]].push_back(i);
		}

		AABB<T> bbox = mesh.getAABB();
		bbox.expand(bbox.extent() * T(0.2));

		if(octree) delete octree;
		octree = new Octree_Hashing<int, Octree_Hashing_AddProcess_TriMesh<Mesh>, T>(bbox, 9);

		Octree_Hashing_AddProcess_TriMesh<Mesh> addProcess(mesh);
		for(unsigned int i=0;i<mesh.numFaces;++i)
			octree->add(int(i), addProcess);

		timer.logTime("ObjectDistanceField3D::init() _ Octree is constructed", log_info);

		set_scale(1.0);
		set_rigid_transformation(arithmetic::Transform<T>());
	}
	const T getPhi(const Vector3<T> &p) const {
		const Vector3<T> pos = transform.getInverse().transformPoint(p);

		Octree_Hashing_MinDistProcess_TriMesh<Mesh> min_dist(mesh);
		octree->lookup(pos, min_dist);

		Vector3<T> normal;
		switch(min_dist.closestPrimitiveType) {
			case Octree_Hashing_MinDistProcess_TriMesh<Mesh>::MINDISTPROCESS_FACE	: normal = mesh.faceNormals[min_dist.closestPrimitiveIdx];	break;
			case Octree_Hashing_MinDistProcess_TriMesh<Mesh>::MINDISTPROCESS_VERTEX	: normal = mesh.normals[min_dist.closestPrimitiveIdx];		break;
			case Octree_Hashing_MinDistProcess_TriMesh<Mesh>::MINDISTPROCESS_EDGE	: normal = mesh.getPseudoNormal_Edge(min_dist.closestPrimitiveIdx, min_dist.closestPrimitiveIdx_edge, t2v); break;

			default : 
				basis::throwError("Error in ObjectDistanceField3D::getPhi()");
		}
		
		Vector3<T> v = pos - min_dist.closestPt;
		return SIGN(dot(v, normal)) * ::abs(v.mag()) * scale_factor;
	}

private : 
	Mesh																	mesh;
	Octree_Hashing<int, Octree_Hashing_AddProcess_TriMesh<Mesh>, T>*		octree;

	Array<vector<unsigned int>>												t2v;

	arithmetic::Transform<T>												rigid_transform;
	T																		scale_factor;

	arithmetic::Transform<T>												transform;
};

template<typename T>
class GridDistanceField3D : public DistanceField3D<T> {
public : 
	typedef typename arithmetic::Grid3<T,T>::index_type index_type;

	GridDistanceField3D() : scale_factor(1.0) {}

	void set_scale(const T s)											 { scale_factor = s;		transform = rigid_transform * arithmetic::scale(scale_factor, scale_factor, scale_factor); }
	void set_rigid_transformation(const arithmetic::Transform<T>& trans) { rigid_transform = trans;	transform = rigid_transform * arithmetic::scale(scale_factor, scale_factor, scale_factor); }

	void set(const arithmetic::Grid3<T,T>& _phi) { phi.resize(_phi); phi = _phi; aabb = phi.getAABB(); }
	void set(const std::string& file, index_type nX, index_type nY, index_type nZ, const std::string& save_filename = ""
		   , bool unitize = true, const arithmetic::Transform<T>& trans = arithmetic::Transform<T>(), bool change_order = false) {

	    Timer timer;

		ObjectDistanceField3D<real> obj_field(file.c_str(), unitize, change_order);
		obj_field.set_rigid_transformation(trans);

		set(obj_field, nX, nY, nZ); 
		if(save_filename.length() != 0)
			save(save_filename.c_str());

		timer.printTime("Time for GridDistanceField3D::set(string& file)");
	}
	void set(const ObjectDistanceField3D<T>& obj, index_type nX, index_type nY, index_type nZ) {
		aabb = obj.getAABB();
		const Vector3<T> ext = aabb.extent();

		phi.init(nX, nY, nZ, ext[0] / (nX-1), ext[1] / (nY-1), ext[2] / (nZ-1));
		phi.set_origin(aabb.min());

#ifdef USE_OPENMP
		int my_end = int(phi.nX);
		int i;
#pragma omp parallel for
#else
		index_type my_end = phi.nX;
		index_type i;
#endif
		for(i=0;i<my_end;++i) 
			for(index_type j=0;j<phi.nY;++j) 
				for(index_type k=0;k<phi.nZ;++k) 
					phi(i,j,k) = obj.getPhi(phi.pos(i,j,k));

		Levelset_Redistancing_PDE_WENO<T> redist(100, 0.2);	
		redist.redistance(phi);

		aabb = phi.getAABB();

		std::cout << "GridDistanceField3D::set(const ObjectDistanceField3D<T>& obj, index_type nX, index_type nY, index_type nZ) Done()" << std::endl;
	}

	const T getPhi(const Vector3<T> &p) const {
		const Vector3<T> pos = transform.getInverse().transformPoint(p);

		if(aabb.isInside(pos))	return phi(pos) * scale_factor;
		else					return 100;
	}

	const bool save(const std::string& filename) const {
		ofstream ofs(StringLibrary::string2wstring(filename.c_str()).c_str(), ios::binary);
		if(!ofs)
			return false;

		phi.exportBinary(ofs);
		return true;
	}

	const bool load(const std::string& filename) {
		ifstream ifs(StringLibrary::string2wstring(filename.c_str()).c_str(), ios::binary);
		if(!ifs)
			return false;

		phi.importBinary(ifs);
		aabb = phi.getAABB();
		set_scale(1.0);
		set_rigid_transformation(arithmetic::Transform<T>());

		return true;
	}


protected : 
	arithmetic::Grid3<T,T>		phi;
	AABB<T>						aabb;

	T							scale_factor;
	arithmetic::Transform<T>	rigid_transform;
	arithmetic::Transform<T>	transform;
};

template<typename T>
class JointDistanceField3D : public DistanceField3D<T> {
public : 
	JointDistanceField3D() {}

	void add_distancefields(DistanceField3D<T>* dists) { fields.push_back(dists); }

	const T getPhi(const Vector3<T> &pos) const {
		T _min = 1.0e+6;
		for(size_t i=0;i<fields.size();++i)
			_min = MIN(_min, fields[i]->getPhi(pos));
		return _min;
	}

protected : 
	std::vector<DistanceField3D<T>*> fields;
};

};		// namespace

#endif