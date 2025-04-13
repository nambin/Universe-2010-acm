#ifndef _GEOMETRYSPACE_TRIMESH_SIMPLE_H_
#define _GEOMETRYSPACE_TRIMESH_SIMPLE_H_

#include <vector>
#include "../../../Basis/Color/RGBA.h"
#include "../../../Basis/Basis/Array.h"
#include "../../Arithmetic/Vector/Vector3.h"
#include "../../Arithmetic/Transform/Transform.h"
#include "../../VirtualSpace/Material/Material.h"
#include "../../Geometry/BiTree/Octree_Hashing.h"
#include "../../Geometry/GeometryLibrary/IntersectionTestLibrary.h"
#include "../../Geometry/GeometryLibrary/DistanceLibrary.h"
#include "OBJStruct.h"
#include "PLYStruct.h"

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

namespace geometryspace {

using std::vector;
using basis::RGBA;
using basis::Array;
using arithmetic::Vector3;
using arithmetic::Transform;
using virtualspace::Material;

template<typename POSVECTOR, typename VERTEXCONTAINER>
class TriMesh_Simple {
public:
	typedef typename POSVECTOR::value_type	real_type;
	typedef typename VERTEXCONTAINER		container_type;

	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	TriMesh_Simple();
	TriMesh_Simple(const TriMesh_Simple &);	
	~TriMesh_Simple();

	void load(const char* const);

	void init();
	void release();

	void set(const TriMesh_Simple &);
	void set(const OBJStruct &obj);
	void set(const PLYStruct &ply);

	void setEpsilon(real_type _eps) { eps = _eps; }

	void exportOBJ(const std::string&, const Transform<real_type>& trans = Transform<real_type>()) const;

	void changeTriangleOrder();

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	void glDraw_FlatSurface(const RGBA&) const;
	void glDraw_SmoothSurface(const RGBA&) const;

	void glDraw_FlatSurface(const Material&) const;
	void glDraw_SmoothSurface(const Material&) const;

	void glDrawLineStrip(const Material&) const;

	//////////////////////////////////////////////////////////////////////////
	// Utility
	void transformMesh(const Transform<real_type> &trans);		// Only works for 3D Mesh
	const real_type getTotalArea() const;
	const real_type getArea(unsigned int face) const;

	const AABB<real_type> getAABB() const;
	const AABB<real_type> getAABB(const arithmetic::Transform<real_type>&) const;

	void unitize();
	void translate(const POSVECTOR &trans);
	void scale(real_type s);	

	//////////////////////////////////////////////////////////////////////////
	// Getter
	const Vector3<real_type>	get_pos(unsigned int vert) const							{ return posVerts[vert];			}
	const Vector3<real_type>	get_pos(unsigned int face, unsigned int idx) const			{ return posVerts[v2f[face*3+idx]]; }

	const unsigned int			get_vert_index(unsigned int face, unsigned int idx) const	{ return v2f[face*3+idx];			}

	//////////////////////////////////////////////////////////////////////////
	// Manipulation
	void makeFaceNormals();
	void makeVertexNormals();
	
	void getVertexNormals(Array<Vector3<real_type>> &) const;

	const POSVECTOR getFaceNormal(unsigned int face);
	const POSVECTOR getPseudoNormal_Vertex(unsigned int iv0, const Array<vector<unsigned int>> &t2n) const;
	const POSVECTOR getPseudoNormal_Edge(unsigned int iv0, unsigned int iv1, const Array<vector<unsigned int>> &t2n) const;

public:
	unsigned int					numVertices;		// # of posVerts
	unsigned int					numFaces;			// # of faces

	VERTEXCONTAINER					posVerts;			// vertex positions
	vector<unsigned int>			v2f;				// three vertex IDs per each triface

	vector<Vector3<real_type>>		faceNormals;		// normal of each triface
	vector<Vector3<real_type>>		normals;

	real_type eps;
};

template <typename POSVECTOR, typename VERTEXCONTAINER> TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::TriMesh_Simple()	{ init();	}
template <typename POSVECTOR, typename VERTEXCONTAINER> TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::TriMesh_Simple(const TriMesh_Simple &mesh) { 
	init();	
	set(mesh);
}
template <typename POSVECTOR, typename VERTEXCONTAINER> TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::~TriMesh_Simple()	{ release();							}
template <typename POSVECTOR, typename VERTEXCONTAINER> void TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::init()			{ release(); eps = (real_type)(1.0e-10);}

template <typename POSVECTOR, typename VERTEXCONTAINER> 
void TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::release() {
	numVertices = 0;
	numFaces = 0;
	posVerts.clear();
	v2f.clear();
	faceNormals.clear();
	normals.clear();
}

template <typename POSVECTOR, typename VERTEXCONTAINER> 
void TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::changeTriangleOrder() {
	for(size_t i=0;i<v2f.size();i+=3) 
		std::swap(v2f[i+1], v2f[i+2]);
}

template <typename POSVECTOR, typename VERTEXCONTAINER> 
void TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::exportOBJ(const std::string& filename, const Transform<real_type>& trans) const {

	std::ofstream ofs(StringLibrary::string2wstring(filename).c_str());

	Vector3<real_type> v;

	ofs << "g default" << endl;
	for(unsigned int i=0;i<numVertices;++i) {
		v = trans.transformPoint(posVerts[i]);
		ofs << "v " << v[0] << " " << v[1] << " " << v[2] << std::endl;
	}

	if(normals.size() != 0) {
		for(unsigned int i=0;i<numVertices;++i) {
			v = trans.transformNormal(normals[i]);
			ofs << "vn " << v[0] << " " << v[1] << " " << v[2] << std::endl;
		}

		for(unsigned int i=0;i<numFaces;++i)
			ofs << "f " << v2f[i*3+0]+1 << "//" << v2f[i*3+0]+1 << " " 
						<< v2f[i*3+1]+1 << "//" << v2f[i*3+1]+1 << " " 
						<< v2f[i*3+2]+1 << "//" << v2f[i*3+2]+1 << std::endl;		
	}
	else {
		for(unsigned int i=0;i<numFaces;++i)
			ofs << "f " << v2f[i*3+0]+1 << " " << v2f[i*3+1]+1 << " " << v2f[i*3+2]+1 << std::endl;
	}	

	ofs.close();
}

template <typename POSVECTOR, typename VERTEXCONTAINER> 
void TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::load(const char * const filename) {

	std::string ext = basis::DataPath::getFilenameExtension(filename);

	if(ext.compare("OBJ")==0 || ext.compare("obj")==0) { OBJStruct obj; obj.loadFromOBJ(filename); set(obj); }	
	else if(ext.compare("PLY")==0 || ext.compare("ply")==0) { PLYStruct ply; ply.loadFromPLY(filename); set(ply); }	
	else
		basis::throwError("Unsupported FileFormat _ TriMesh_Simple::load(..)\n");

	return;
}

template <typename POSVECTOR, typename VERTEXCONTAINER>
void TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::set(const TriMesh_Simple &mesh) {
	release();

	numVertices = mesh.numVertices;
	numFaces = mesh.numFaces;

	posVerts.resize(mesh.posVerts.size());
	v2f.resize(mesh.v2f.size());
	faceNormals.resize(mesh.faceNormals.size());
	normals.resize(mesh.normals.size());

	for( size_t i=0; i<mesh.posVerts.size(); ++i ) 		posVerts[i] = mesh.posVerts[i];
	for( size_t i=0; i<mesh.v2f.size(); ++i )			v2f[i] = mesh.v2f[i];
	for( size_t i=0; i<mesh.faceNormals.size(); ++i ) 	faceNormals[i] =  mesh.faceNormals[i];
	for( size_t i=0; i<mesh.normals.size(); ++i )		normals[i] = mesh.normals[i];
}

template <typename POSVECTOR, typename VERTEXCONTAINER>
void TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::set(const OBJStruct &obj) {
	release();

	numVertices = obj.numVertices;
	numFaces = obj.numTriangles;

	posVerts.resize(numVertices);
	v2f.resize(numFaces * 3);

	for( size_t i=0; i<posVerts.size(); ++i ) 						posVerts[i] = obj.posVerts[i];
	for( size_t i=0; i<numFaces; ++i ) for( size_t j=0;j<3;++j)		v2f[i*3+j] = obj.triangles[i].vIndex[j];
	
	makeFaceNormals();
	makeVertexNormals();
}

template <typename POSVECTOR, typename VERTEXCONTAINER>
void TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::set(const PLYStruct &ply) {
	release();

	numVertices = ply.numVertices;
	numFaces = ply.numTriangles;

	v2f.resize(numFaces*3);			for(unsigned int i=0;i<numFaces;++i) for(unsigned int j=0;j<3;++j) v2f[i*3+j] = ply.triangles[i].vIndex[j];
	posVerts.resize(numVertices);	for(unsigned int i=0;i<numVertices;++i) posVerts[i] = ply.vertices[i].pos;

	makeFaceNormals();
	makeVertexNormals();
}

template <typename POSVECTOR, typename VERTEXCONTAINER>
void TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::glDraw_FlatSurface(const RGBA& color) const {
	Material mat;
	mat.setAmbient(color);
	mat.setDiffuse(color);
	
	glDraw_FlatSurface(mat);
}

template <typename POSVECTOR, typename VERTEXCONTAINER>
void TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::glDraw_SmoothSurface(const RGBA& color) const {
	Material mat;
	mat.setAmbient(color);
	mat.setDiffuse(color);

	glDraw_SmoothSurface(mat);
}

template <typename POSVECTOR, typename VERTEXCONTAINER>
void TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::glDraw_FlatSurface(const Material& mat) const {

	glEnable(GL_LIGHTING);
	glShadeModel(GL_FLAT);

	glDisable(GL_COLOR_MATERIAL);
	mat.glBindMaterial();

	glBegin( GL_TRIANGLES );

	for(unsigned int i=0;i<numFaces;++i)	{
		unsigned int i0 = v2f[i*3  ];
		unsigned int i1 = v2f[i*3+1];
		unsigned int i2 = v2f[i*3+2];

		const Vector3<real_type> *V0, *V1, *V2;
		V0 = &posVerts[v2f[i*3]  ];
		V1 = &posVerts[v2f[i*3+1]];
		V2 = &posVerts[v2f[i*3+2]];

		glNormal3f(faceNormals[i][0], faceNormals[i][1], faceNormals[i][2]);

		glVertex3f((*V0)[0], (*V0)[1], (*V0)[2]);
		glVertex3f((*V1)[0], (*V1)[1], (*V1)[2]);
		glVertex3f((*V2)[0], (*V2)[1], (*V2)[2]);
	}

	glEnd();
}

template <typename POSVECTOR, typename VERTEXCONTAINER>
void TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::glDraw_SmoothSurface(const Material& mat) const {
	if(normals.empty()) {
		glDraw_FlatSurface(mat);
		return;
	}

	glEnable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);

	glDisable(GL_COLOR_MATERIAL);
	mat.glBindMaterial();	

	glBegin( GL_TRIANGLES );

	for(unsigned int i=0;i<numFaces;++i)	{
		unsigned int i0 = v2f[i*3  ];
		unsigned int i1 = v2f[i*3+1];
		unsigned int i2 = v2f[i*3+2];

		const Vector3<real_type> *V0, *V1, *V2, *N0, *N1, *N2;
		V0 = &posVerts[v2f[i*3]  ];
		V1 = &posVerts[v2f[i*3+1]];
		V2 = &posVerts[v2f[i*3+2]];
		N0 = &normals[v2f[i*3]  ];
		N1 = &normals[v2f[i*3+1]];
		N2 = &normals[v2f[i*3+2]];

		glNormal3f((*N0)[0], (*N0)[1], (*N0)[2]);
		glVertex3f((*V0)[0], (*V0)[1], (*V0)[2]);
		glNormal3f((*N1)[0], (*N1)[1], (*N1)[2]);
		glVertex3f((*V1)[0], (*V1)[1], (*V1)[2]);
		glNormal3f((*N2)[0], (*N2)[1], (*N2)[2]);
		glVertex3f((*V2)[0], (*V2)[1], (*V2)[2]);
	}

	glEnd();
}

template <typename POSVECTOR, typename VERTEXCONTAINER>
void TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::glDrawLineStrip(const Material& mat) const {
//	glEnable(GL_LIGHTING);
//	glShadeModel(GL_SMOOTH);

//	glDisable(GL_COLOR_MATERIAL);
//	mat.glBindMaterial();	

	glDisable(GL_LIGHTING);
	glColor4f(mat.getDiffuse()[0], mat.getDiffuse()[1], mat.getDiffuse()[2], mat.getDiffuse()[3]);

	glLineWidth(0.5f);
	glBegin( GL_TRIANGLE_STRIP );

	for(unsigned int i=0;i<numFaces;++i)	{
		unsigned int i0 = v2f[i*3  ];
		unsigned int i1 = v2f[i*3+1];
		unsigned int i2 = v2f[i*3+2];

		const Vector3<real_type> *V0, *V1, *V2;
		V0 = &posVerts[v2f[i*3]  ];
		V1 = &posVerts[v2f[i*3+1]];
		V2 = &posVerts[v2f[i*3+2]];

		glVertex3f((*V0)[0], (*V0)[1], (*V0)[2]);
		glVertex3f((*V1)[0], (*V1)[1], (*V1)[2]);
		glVertex3f((*V2)[0], (*V2)[1], (*V2)[2]);
	}

	glEnd();
}

template <typename POSVECTOR, typename VERTEXCONTAINER>
void TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::transformMesh(const arithmetic::Transform<real_type> &trans) {
	for(size_t i=0;i<posVerts.size();++i)
		posVerts[i] = trans.transformPoint(posVerts[i]);

	if(faceNormals.size() != 0) {
		for(size_t i=0;i<faceNormals.size();++i)
			faceNormals[i] = trans.transformNormal(faceNormals[i]);
	}
	if(normals.size() != 0) {
		for(size_t i=0;i<normals.size();++i)
			normals[i] = trans.transformNormal(normals[i]);
	}
}

template <typename POSVECTOR, typename VERTEXCONTAINER>
const typename TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::real_type TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::getTotalArea() const {
	real_type area = 0;
	for(unsigned int i=0;i<numFaces;++i)
		area += getArea(i);
	return area;
}

template <typename POSVECTOR, typename VERTEXCONTAINER>
const typename TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::real_type TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::getArea(unsigned int triface) const {
	Vector3<real_type> V01 = posVerts[ v2f[3*triface+1] ] - posVerts[ v2f[3*triface  ] ];
	Vector3<real_type> V02 = posVerts[ v2f[3*triface+2] ] - posVerts[ v2f[3*triface  ] ];

	return arithmetic::mag(cross(V01, V02)) * 0.5f;
}

template <typename POSVECTOR, typename VERTEXCONTAINER>
const AABB<typename POSVECTOR::value_type> TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::getAABB() const {
	AABB<typename POSVECTOR::value_type> aabb;

	for(unsigned int i=0;i<posVerts.size();++i)
		aabb.surround(posVerts[i]);

	return aabb;
}

template <typename POSVECTOR, typename VERTEXCONTAINER>
const AABB<typename POSVECTOR::value_type> TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::getAABB(const arithmetic::Transform<real_type>& trans) const {
	AABB<typename POSVECTOR::value_type> aabb;

	for(unsigned int i=0;i<posVerts.size();++i)
		aabb.surround(trans.transformPoint(posVerts[i]));

	return aabb;
}

template <typename POSVECTOR, typename VERTEXCONTAINER>
void TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::unitize() {
	real_type maxx, minx, maxy, miny, maxz, minz;

	// get the max/min
	maxx = minx = posVerts[0][0];
	maxy = miny = posVerts[0][1];
	maxz = minz = posVerts[0][2];

	for(unsigned int i = 1; i < numVertices; ++i) {
		if (maxx < posVerts[i][0])	maxx = posVerts[i][0];
		if (minx > posVerts[i][0])	minx = posVerts[i][0];

		if (maxy < posVerts[i][1])	maxy = posVerts[i][1];
		if (miny > posVerts[i][1])	miny = posVerts[i][1];

		if (maxz < posVerts[i][2])	maxz = posVerts[i][2];
		if (minz > posVerts[i][2])	minz = posVerts[i][2];
	}

	// calculate model width, height, and depth 
	real_type w = abs(maxx-minx);
	real_type h = abs(maxy-miny);
	real_type d = abs(maxz-minz);

	// calculate center of the model
	real_type cx = (maxx + minx) * 0.5f;
	real_type cy = (maxy + miny) * 0.5f;
	real_type cz = (maxz + minz) * 0.5f;

	// calculate unitizing scale factor 
	real_type s = 1.0f / basis::MAX3(w, h, d);

	// translate around center then scale 
	translate(-POSVECTOR(cx, cy, cz));
	scale(s);
}

template <typename POSVECTOR, typename VERTEXCONTAINER>
void TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::translate(const POSVECTOR &trans) {
	for(unsigned int i = 0; i < numVertices; ++i) 
		posVerts[i] += trans;
}

template <typename POSVECTOR, typename VERTEXCONTAINER>
void TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::scale(real_type s)	{
	for(unsigned int i = 0; i < numVertices; ++i) {
		posVerts[i][0] *= s;
		posVerts[i][1] *= s;
		posVerts[i][2] *= s;
	}
}

template <typename POSVECTOR, typename VERTEXCONTAINER>
const POSVECTOR TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::getFaceNormal(unsigned int triface) {

	Vector3<real_type> V01 = posVerts[ v2f[3*triface+1] ] - posVerts[ v2f[3*triface  ] ];
	Vector3<real_type> V02 = posVerts[ v2f[3*triface+2] ] - posVerts[ v2f[3*triface  ] ];

	Vector3<real_type> normal = cross(V01, V02);
	real_type mag = arithmetic::mag(normal);

	if(IsAlmostZero(mag, eps))	{
		std::cerr << "Triface_Simple Normal Magnitude : " << mag << std::endl;
		std::cerr << "Normal : " << normal << std::endl;
		std::cerr << "V01 : " << V01 << std::endl;
		std::cerr << "V02 : " << V02 << std::endl;
	}
	else 
		normal /= mag;

	return normal;
}

template <typename POSVECTOR, typename VERTEXCONTAINER>
void TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::makeFaceNormals()	{

	faceNormals.clear();
	faceNormals.resize(numFaces);

	for( unsigned int i=0; i<numFaces; ++i)	
		faceNormals[i] = getFaceNormal(i);
}

template <typename POSVECTOR, typename VERTEXCONTAINER>
void TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::makeVertexNormals()	{

	Array<Vector3<real_type>> vertNormals(posVerts.size());
	getVertexNormals(vertNormals);

	normals.resize(posVerts.size());
	for(size_t i=0;i<vertNormals.size();++i)
		normals[i] = vertNormals[i];
}

template <typename POSVECTOR, typename VERTEXCONTAINER>
void TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::getVertexNormals(Array<Vector3<real_type>> &_vertNormals) const {
	assert(_vertNormals.size() == posVerts.size());

	if(faceNormals.size() != numFaces)	{
		throwError("Triface_Simple faceNormals are not updated");
		return;
	}

	Array<vector<unsigned int>> t2v(numVertices);
	for(unsigned int i=0;i<numFaces;i++)	{
		t2v[v2f[i*3  ]].push_back(i);
		t2v[v2f[i*3+1]].push_back(i);
		t2v[v2f[i*3+2]].push_back(i);
	}

	Vector3<real_type> normal;
	for(unsigned int i=0;i<numVertices;i++)	{

		normal = getPseudoNormal_Vertex(i, t2v);

		real_type mag = arithmetic::mag(normal);
		if(IsAlmostZero(mag, eps))	{
//			std::cout << "In PseudoNormal on Vertex Magnitude : " << mag << std::endl;
//			std::cout << "Vertex Idx : " << i << " _ Normal : " << normal << std::endl;			
		}
		else
			normal /= mag;

		_vertNormals[i] = normal;
	}
}

template <typename POSVECTOR, typename VERTEXCONTAINER>
const POSVECTOR TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::getPseudoNormal_Vertex(unsigned int iv0, const Array<vector<unsigned int>> &t2n) const {

	POSVECTOR normal(0.0f);

	const Vector3<real_type> *v0;
	Vector3<real_type> v[2];

	v0 = &posVerts[iv0];

	for(size_t i=0;i<t2n[iv0].size();i++)	{

		unsigned int idxTemp = 0;
		unsigned int face = t2n[iv0][i];

		for(unsigned int j=0;j<3;j++)	{	
			unsigned int iv = v2f[face*3+j];
			if( iv != iv0 )	
				v[idxTemp++] = posVerts[iv];
		}

		v[0] -= (*v0);
		v[1] -= (*v0);
		real_type _angle = angle(v[0], v[1]);

		normal += (faceNormals[face] * _angle);

	}

	real_type mag_ = mag(normal);
	if(IsAlmostZero(mag_, eps))	{
//		std::cout << "PseudoNormal On Vertex Magnitude : " << mag_ << std::endl;
//		std::cout << "Normal : " << normal << std::endl;
	}

	return normal;
}

template <typename POSVECTOR, typename VERTEXCONTAINER>
const POSVECTOR TriMesh_Simple<POSVECTOR,VERTEXCONTAINER>::getPseudoNormal_Edge(unsigned int iv0, unsigned int iv1, const Array<vector<unsigned int>> &t2n) const {

	unsigned int nTriface = 0;
	unsigned int iTriface[10];
	for(size_t i=0;i<t2n[iv0].size();i++)	{
		for(unsigned int j=0;j<3;j++)	{
			if(v2f[t2n[iv0][i]*3+j] == iv1)	{
				iTriface[nTriface++] = t2n[iv0][i];
			}
		}
	}

	if(nTriface > 2)	{
//		std::cout << "In get PseudoNormal_Edge" << std::endl;
//		std::cout << iv0 << " " << iv1 << " " << nTriface << std::endl;
		for(unsigned int i=0;i<nTriface;i++)	{
			std::cout << iTriface[i] << std::endl;

			unsigned int iVertex[3];

			iVertex[0] = v2f[3*iTriface[i]  ];
			iVertex[1] = v2f[3*iTriface[i]+1];
			iVertex[2] = v2f[3*iTriface[i]+2];

			Vector3<real_type> V0 = posVerts[iVertex[0]];
			Vector3<real_type> V1 = posVerts[iVertex[1]];
			Vector3<real_type> V2 = posVerts[iVertex[2]];

//			std::cout << iVertex[0] << " : " << V0 << std::endl;
//			std::cout << iVertex[1] << " : " << V1 << std::endl;
//			std::cout << iVertex[2] << " : " << V2 << std::endl;
		}
	}

	POSVECTOR normal(0.0f);
	for(unsigned int i=0;i<nTriface;i++)	
		normal += faceNormals[iTriface[i]];

	real_type mag_ = mag(normal);
	if(IsAlmostZero(mag_, eps))	{
		//std::cout << "PseudoNormal On Edge Magnitude : " << mag_ << std::endl;
		//std::cout << "Number of Triface : " << nTriface << std::endl;
		//std::cout << "Vertex Idx : " << iv0 << " " << iv1 << std::endl;
		//std::cout << "Normal : " << normal << std::endl;
	}

	return normal;
}

////////////////////////////////////////
////////////////////////////////////////
// Functor for Octree_Hashing  /////////
////////////////////////////////////////
////////////////////////////////////////

template<class TRIMESH>
class Octree_Hashing_AddProcess_TriMesh {
public :
	typedef typename TRIMESH::real_type T;

	Octree_Hashing_AddProcess_TriMesh(const TRIMESH& _mesh) : mesh(_mesh) {}

	OCTREE_ADD_OPTION operator()(int idx, const AABB<T> &bbox) const {		
		if(IntersectionTestLibrary::testTriangleAABB(mesh.get_pos(idx, 0), mesh.get_pos(idx, 1), mesh.get_pos(idx, 2), bbox.min(), bbox.max()))	return OCTREE_REFINE;
		else																																	return OCTREE_NONE;
	}

public :
	const TRIMESH& mesh;
};

template<class TRIMESH>
class Octree_Hashing_MinDistProcess_TriMesh {
public : 
	enum { MINDISTPROCESS_VERTEX = 0, MINDISTPROCESS_EDGE, MINDISTPROCESS_FACE };

	typedef typename TRIMESH::real_type T;

public : 
	Octree_Hashing_MinDistProcess_TriMesh(const TRIMESH& _mesh) : mesh(_mesh) { init(); }
	void init() { minSqDist = 1.0e+8; closestTriangleIdx = -1; }

	const bool operator()(const Vector3<T> &p, const AABB<T> &bbox) const { return (minSqDist > DistanceLibrary::SqDistPointAABB(p, bbox.min(), bbox.max())); }

	void operator()(const Vector3<T> &p, int idx) const {
		Vector3<T> _closestPt; DistanceLibrary::Closest_Point_Type _closestPtType;

		const Vector3<T> _bary = DistanceLibrary::ClosestPtPointTriangle(_closestPt, _closestPtType, p, mesh.get_pos(idx, 0), mesh.get_pos(idx, 1), mesh.get_pos(idx, 2));
		const T _sqDist = sqrlength(_closestPt, p);

		if(minSqDist > _sqDist) {
			closestPt = _closestPt;
			closestTriangleIdx = idx;
			minSqDist = _sqDist;
			baryCentric = _bary;

			switch(_closestPtType) {
				case DistanceLibrary::CLOSEST_PT_EDGE0 : 
				case DistanceLibrary::CLOSEST_PT_EDGE1 : 
				case DistanceLibrary::CLOSEST_PT_EDGE2 : 					
				{
					closestPrimitiveType = MINDISTPROCESS_EDGE;

					int lc_idx = _closestPtType - DistanceLibrary::CLOSEST_PT_EDGE0;

					closestPrimitiveIdx			= mesh.v2f[3*idx + (lc_idx + 1) % 3];
					closestPrimitiveIdx_edge	= mesh.v2f[3*idx + (lc_idx + 2) % 3];

					break;
				}

				case DistanceLibrary::CLOSEST_PT_NODE0 : 
				case DistanceLibrary::CLOSEST_PT_NODE1 : 
				case DistanceLibrary::CLOSEST_PT_NODE2 : 
				{
					closestPrimitiveType = MINDISTPROCESS_VERTEX;

					int lc_idx = _closestPtType - DistanceLibrary::CLOSEST_PT_NODE0;

					closestPrimitiveIdx = mesh.v2f[3*idx + lc_idx];
					closestPrimitiveIdx_edge = -1;

					break;
				}
				default	: 

					closestPrimitiveType = MINDISTPROCESS_FACE;

					closestPrimitiveIdx = idx;
					closestPrimitiveIdx_edge = -1;

					break;
			}

		}		
	}
	const T getMinSqDist() const						{ return minSqDist;				}
	const T getMinDist() const							{ return sqrt(minSqDist);		}

	const Vector3<T> getClosestPoint() const			{ return closestPt;				}
	const Vector3<T> getBaryCentricCoordinate() const	{ return baryCentric;			}

	const int getClosestTriangleIdx() const				{ return closestTriangleIdx;	}
	const int getClosestPrimitiveType() const			{ return closestPrimitiveType;	}
	const int getClosestPrimitiveIdx() const			{ return closestPrimitiveIdx;	}

public : 
	const TRIMESH&			mesh;
	mutable T				minSqDist;
	mutable Vector3<T>		closestPt;
	mutable int				closestTriangleIdx;

	mutable Vector3<T>		baryCentric;

	mutable int				closestPrimitiveType;
	mutable int				closestPrimitiveIdx;
	mutable int				closestPrimitiveIdx_edge;
};

};		// namespace geometryspace

#endif