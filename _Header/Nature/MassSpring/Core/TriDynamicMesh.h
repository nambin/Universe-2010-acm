#ifndef _MASS_SPRING_TRI_DYNAMICMESH_H_
#define _MASS_SPRING_TRI_DYNAMICMESH_H_

#include "../../../Basis/Color/RGBA.h"
#include "../../../Principle/Arithmetic/Transform/Transform.h"
#include "../../../Principle/Geometry/Surface/TriMesh_Simple.h"
#include "../../../Principle/Geometry/Surface/TriMesh.h"
#include "../Core/ClothLibrary.h"
#include "../Core/Constraints.h"
#include "../Core/MassSpringSystemMatrix.h"
#include "../Core/DynamicMesh.h"

namespace springspace {

using basis::RGBA;
using basis::Array2;
using arithmetic::Vector2;
using geometryspace::TriMesh;
using geometryspace::TriMesh_Simple;

template<typename T>
class TriDynamicMesh : public DynamicMesh, public TriMesh_Simple<Vector3<T>, Vec3Array<T>> {
public :
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	TriDynamicMesh();
	template<typename CONTAINER> TriDynamicMesh(const TriMesh_Simple<Vector2<T>, CONTAINER>&);
	TriDynamicMesh(const TriMesh<Vector2<T>>&);
	~TriDynamicMesh();

	void postProcessing();

	void init();
	void initpos();
	template<typename CONTAINER> void initpos(const TriMesh_Simple<Vector2<T>, CONTAINER>&);
	void initpos(const TriMesh<Vector2<T>>&);

	void make_WidthAndHeight();

	void setDensity(T);
	void setMass(T);
	void setParticleMass(T);	

	float width() const  { return float(Lx); }
	float length() const { return float(Ly); }
	index_type size() const  { return numVertices; }

	//////////////////////////////////////////////////////////////////////////
	// Spring 
	void make_StretchSpring(SpringFunctor<T> *);
	void make_ShearSpring(SpringFunctor<T> *);
	void make_BendingSpring(SpringFunctor<T> *);

	void allocate_SystemMatrix(MassSpringSystemMatrix<T>&) const;

	//////////////////////////////////////////////////////////////////////////
	// Constraint
	void fixVertex(size_t i);

	//////////////////////////////////////////////////////////////////////////
	// Utility
	void makeV2V();
	void transformMesh(const arithmetic::Transform<T> &trans);
	const Vector3<T> centerOfMass() const;

	void calcDensity();
	const T getDensity() const;
	const T maxStrain() const;

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	void glDraw_FlatSurface(const Material& mf, const Material& mb) const;
	void glDraw_SmoothSurface(const Material& mf, const Material& mb) const;

	void glDraw_Particle(const RGBA& c0, const RGBA& cConstraint, GLfloat size) const;
	void glDraw_StretchSpring(const RGBA&, GLfloat width) const;
	void glDraw_ShearSpring(const RGBA&, GLfloat width) const;
	void glDraw_BendingSpring(const RGBA&, GLfloat width) const;

public :
	T							Lx, Ly;

	Array<Vector2<T>>			pos_rest;

	Vec3Array<T>&				pos;
	Vec3Array<T>				pos0;
	Vec3Array<T>				vel, vel0;

	Array<T>					mass;
	T							density;

	SpringFunctor<T>*			stretchFunc;
	SpringFunctor<T>*			shearFunc; 
	SpringFunctor<T>*			bendingFunc;
	Array<Spring<T>>			stretchSpring;		
	Array<Spring<T>>			shearSpring;	
	Array<Spring<T>>			bendingSpring;	

	Constraints<T>				constraint;

	Array<vector<Spring<T>*>>	s2v;
	Array<vector<size_t>>		v2v;
};

template<typename T>
TriDynamicMesh<T>::TriDynamicMesh() 
: TriMesh_Simple<Vector3<T>, Vec3Array<T>>(), pos(posVerts), stretchFunc(NULL), shearFunc(NULL), bendingFunc(NULL)
{}

template<typename T> template<typename CONTAINER> 
TriDynamicMesh<T>::TriDynamicMesh(const geometryspace::TriMesh_Simple<Vector2<T>, CONTAINER>& mesh) 
: TriMesh_Simple<Vector3<T>, Vec3Array<T>>(), pos(posVerts), stretchFunc(NULL), shearFunc(NULL), bendingFunc(NULL)
{ initpos(mesh); }

template<typename T> 
TriDynamicMesh<T>::TriDynamicMesh(const geometryspace::TriMesh<Vector2<T>>& mesh) 
: TriMesh_Simple<Vector3<T>, Vec3Array<T>>(), pos(posVerts), stretchFunc(NULL), shearFunc(NULL), bendingFunc(NULL)
{ initpos(mesh); }

template<typename T> 
TriDynamicMesh<T>::~TriDynamicMesh() {
	if(stretchFunc)	delete stretchFunc;	stretchFunc = NULL;
	if(shearFunc)	delete shearFunc;	shearFunc = NULL;
	if(bendingFunc)	delete bendingFunc;	bendingFunc = NULL;
}

template<typename T>
void TriDynamicMesh<T>::postProcessing() {
	makeFaceNormals();
	makeVertexNormals();
}

template<typename T>
void TriDynamicMesh<T>::init() {
	TriMesh_Simple<Vector3<T>, Vec3Array<T>>::init();

	pos_rest.clear(); 
	pos.clear(); pos0.clear();
	vel.clear(); vel0.clear();	
	mass.clear();

	stretchSpring.clear(); shearSpring.clear(); bendingSpring.clear();

	if(stretchFunc)	delete stretchFunc;	stretchFunc = NULL;
	if(shearFunc)	delete shearFunc;	shearFunc = NULL;
	if(bendingFunc)	delete bendingFunc;	bendingFunc = NULL;

	constraint.clear();

	s2v.clear();
	v2v.clear();
}

template<typename T>
void TriDynamicMesh<T>::initpos() {
	for(size_t i=0;i<pos_rest.size();++i) 
		pos[i] = pos0[i] = Vector3<T>(pos_rest[i][0], pos_rest[i][1], 0);

	makeFaceNormals();
	makeVertexNormals();
}

template<typename T> template<typename CONTAINER> 
void TriDynamicMesh<T>::initpos(const TriMesh_Simple<Vector2<T>, CONTAINER>& mesh) {
	init();

	numVertices = mesh.numVertices;
	numFaces = mesh.numFaces;

	pos.resize(numVertices);
	pos0.resize(numVertices);
	vel.resize(numVertices);
	vel0.resize(numVertices);
	pos_rest.resize(numVertices);
	v2f.resize(numFaces*3);

	mass.resize(numVertices);	mass.set(1.0f);

	constraint.init(numVertices);	

	s2v.resize(numVertices);

	for(size_t i=0;i<mesh.vertices.size();++i) 	pos_rest[i] = mesh.vertices[i];
	for(size_t i=0;i<mesh.v2f.size();++i) 	v2f[i] = mesh.v2f[i];
	initpos();	

	calcDensity();
	makeV2V();
	make_WidthAndHeight();	
}

template<typename T>
void TriDynamicMesh<T>::make_WidthAndHeight() {
	T minX = 1.0e+6, minY = 1.0e+6, maxX = -1.0e+6, maxY = -1.0e+6;
	for(size_t i=0;i<pos_rest.size();++i) {
		if(minX > pos_rest[i][0]) minX = pos_rest[i][0];
		if(maxX < pos_rest[i][0]) maxX = pos_rest[i][0];
		if(minY > pos_rest[i][1]) minY = pos_rest[i][1];
		if(maxY < pos_rest[i][1]) maxY = pos_rest[i][1];
	}
	Lx = maxX - minX;
	Ly = maxY - minY;
}

template<typename T>
void TriDynamicMesh<T>::initpos(const TriMesh<Vector2<T>>& mesh) {
	init();

	numVertices = mesh.numVertices;
	numFaces = mesh.numFaces;

	pos.resize(numVertices);
	pos0.resize(numVertices);
	vel.resize(numVertices);
	vel0.resize(numVertices);
	pos_rest.resize(numVertices);
	v2f.resize(numFaces*3);

	mass.resize(numVertices);	mass.set(1.0f);

	constraint.init(numVertices);

	s2v.resize(numVertices);

	for(size_t i=0;i<mesh.posVerts.size();++i) 	pos_rest[i] = mesh.posVerts[i];
	for(size_t i=0;i<mesh.v2f.size();++i) 		v2f[i] = mesh.v2f[i];
	initpos();

	calcDensity();
	makeV2V();
	make_WidthAndHeight();
}

template<typename T>
void TriDynamicMesh<T>::setMass(T m) {
	m /= mass.size();
	setParticleMass(m);
}

template<typename T>
void TriDynamicMesh<T>::setParticleMass(T m) {
	for(size_t i=0;i<mass.size();++i)
		mass[i] = m;
	calcDensity();
}

template<typename T>
void TriDynamicMesh<T>::setDensity(T den) {
	T area = getTotalArea();
	T massSum = area * den;
	setMass(massSum);
}

template<typename T>
void TriDynamicMesh<T>::make_StretchSpring(SpringFunctor<T> *func) {
	if(!stretchSpring.empty()) {
		std::cout << "Stretch Springs are already made" << std::endl;
		return;
	}

	stretchFunc = func;

	TriMesh<Vector3<T>> mesh;
	mesh.set(*this);
	mesh.makeConnectivityData();

	std::vector<Spring<T>> springs;
	for(size_t i=0;i<mesh.v2e.dimX();++i) 
		springs.push_back(Spring<T>(mesh.v2e(i,0), mesh.v2e(i,1), arithmetic::length(pos[mesh.v2e(i,0)], pos[mesh.v2e(i,1)])));

	stretchSpring.resize(springs.size());
	stretchSpring.set(springs);	
	
	for(size_t i=0;i<stretchSpring.size();++i) {
		s2v[stretchSpring[i].idx0].push_back(&stretchSpring[i]);
		s2v[stretchSpring[i].idx1].push_back(&stretchSpring[i]);
	}
}

template<typename T>
void TriDynamicMesh<T>::make_ShearSpring(SpringFunctor<T> *func) {
}

template<typename T>
void TriDynamicMesh<T>::make_BendingSpring(SpringFunctor<T> *func) {
}

template<typename T>
void TriDynamicMesh<T>::allocate_SystemMatrix(MassSpringSystemMatrix<T>& mat) const {
	mat.clear();

	mat.n = size();
	mat.diag.resize(mat.n);

	size_t num = 0;
	for(size_t i=0;i<s2v.size();++i)
		num += s2v[i].size();

	mat.value.resize(num);
	mat.colindex.resize(num);
	mat.rowstart.resize(mat.n+1);

	mat.rowstart[0] = 0;
	for(size_t i=0;i<s2v.size();++i) {
		mat.rowstart[i+1] = mat.rowstart[i] + s2v[i].size();
		for(size_t j=0;j<s2v[i].size();++j) {
			if(i == s2v[i][j]->idx0) { mat.colindex[mat.rowstart[i]+j] = s2v[i][j]->idx1; }
			else					 { mat.colindex[mat.rowstart[i]+j] = s2v[i][j]->idx0; }
		}
	}
}

template<typename T>
void TriDynamicMesh<T>::fixVertex(size_t _idx) {
	constraint.addConstraint(_idx, Vector3<T>(1,0,0));
	constraint.addConstraint(_idx, Vector3<T>(0,1,0));
	constraint.addConstraint(_idx, Vector3<T>(0,0,1));
	constraint.setVelocity(_idx, Vector3<T>(0));
	constraint.setHardConstraint(_idx);
}

template<typename T>
void TriDynamicMesh<T>::makeV2V() {
	v2v.clear();

	TriMesh<Vector3<T>> mesh;
	mesh.set(*this);
	mesh.makeConnectivityData();

	v2v.resize(numVertices);
	//for(size_t i=0;i<numVertices;++i) {
	//	v2v[i].resize(mesh.v2v[i].size());
	//	for(size_t k=0;k<mesh.v2v[i].size();++k) 
	//		v2v[i][k] = mesh.v2v[i][k];
	//}

	for(int i=0;i<mesh.numEdges;++i) {
		v2v[mesh.v2e[i*2+1]].push_back(mesh.v2e[i*2]);
		v2v[mesh.v2e[i*2]].push_back(mesh.v2e[i*2+1]);
	}

//	cout << v2v[0].size() << endl;
//	cout << mesh.v2v[0].size() << endl;
//	cout << mesh.e2v[0].size() << endl;
//	cout << mesh.f2v[0].size() << endl;
}

template<typename T>
void TriDynamicMesh<T>::transformMesh(const arithmetic::Transform<T> &trans) {
	TriMesh_Simple<Vector3<T>, Vec3Array<T>>::transformMesh(trans);
	for(size_t i=0;i<pos.size();++i)
		pos0(i) = pos(i);

	calcDensity();
}

template<typename T>
const Vector3<T> TriDynamicMesh<T>::centerOfMass() const {
	Vector3<T> com(0);
	T massSum = 0;
	for(size_t i=0;i<pos.size();++i) {
		com += pos(i) * mass(i);
		massSum += mass(i);
	}
	com /= massSum;
	return com;
}

template<typename T>
void TriDynamicMesh<T>::calcDensity() {
	T massSum = 0;
	for(size_t i=0;i<mass.size();++i)
		massSum += mass[i];

	T area = getTotalArea();
	density = massSum / area;

	cout << massSum << " " << area << endl;
}

template<typename T>
const T TriDynamicMesh<T>::getDensity() const {
	return density;
}

template<typename T>
const T TriDynamicMesh<T>::maxStrain() const {
	return 1.0f;
}

template<typename T>
void TriDynamicMesh<T>::glDraw_FlatSurface(const Material& mf, const Material& mb) const {	
	glShadeModel(GL_FLAT);

	glEnable(GL_LIGHTING);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);

	glDisable(GL_CULL_FACE);
	glDisable(GL_COLOR_MATERIAL);	

	mf.glBind_FrontMaterial();
	mb.glBind_BackMaterial();

	glBegin( GL_TRIANGLES );

	const Vector3<T> *V0, *V1, *V2;
	Vector3<T> normal;
	for(size_t i=0;i<numFaces;++i)	{
		V0 = &pos(v2f[i*3]);	V1 = &pos(v2f[i*3+1]);	V2 = &pos(v2f[i*3+2]);
		normal = unitVector(cross(*V1-*V0, *V2-*V0));

		glNormal3f(normal[0], normal[1], normal[2]);

		glVertex3f((*V0)[0], (*V0)[1], (*V0)[2]);
		glVertex3f((*V1)[0], (*V1)[1], (*V1)[2]);
		glVertex3f((*V2)[0], (*V2)[1], (*V2)[2]);
	}

	glEnd();

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 0);
}

template<typename T>
void TriDynamicMesh<T>::glDraw_SmoothSurface(const Material& mf, const Material& mb) const {
	glShadeModel(GL_SMOOTH);

	glEnable(GL_LIGHTING);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);

	glDisable(GL_CULL_FACE);
	glDisable(GL_COLOR_MATERIAL);	

	mf.glBind_FrontMaterial();
	mb.glBind_BackMaterial();

	glBegin( GL_TRIANGLES );

	const Vector3<T> *V0, *V1, *V2;
	for(size_t i=0;i<numFaces;++i)	{
		V0 = &pos(v2f[i*3]);	V1 = &pos(v2f[i*3+1]);	V2 = &pos(v2f[i*3+2]);

		glNormal3f(normals[v2f[i*3]][0], normals[v2f[i*3]][1], normals[v2f[i*3]][2]);
		glVertex3f((*V0)[0], (*V0)[1], (*V0)[2]);
		glNormal3f(normals[v2f[i*3+1]][0], normals[v2f[i*3+1]][1], normals[v2f[i*3+1]][2]);
		glVertex3f((*V1)[0], (*V1)[1], (*V1)[2]);
		glNormal3f(normals[v2f[i*3+2]][0], normals[v2f[i*3+2]][1], normals[v2f[i*3+2]][2]);
		glVertex3f((*V2)[0], (*V2)[1], (*V2)[2]);
	}

	glEnd();

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 0);
}

template<typename T>
void TriDynamicMesh<T>::glDraw_Particle(const RGBA& c0, const RGBA& c1, GLfloat size) const {
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_LIGHTING);

	GLUquadricObj *qobj = gluNewQuadric();

	gluQuadricDrawStyle(qobj, (GLenum) GLU_FILL);
	gluQuadricNormals(qobj, (GLenum) GLU_SMOOTH);

	for(size_t i=0;i<pos.size();++i) {
		if(constraint.isConstrained(i))	glColor4f(c1[0], c1[1], c1[2], c1[3]);
		else							glColor4f(c0[0], c0[1], c0[2], c0[3]);

		glTranslatef(pos[i][0], pos[i][1], pos[i][2]);		
		gluSphere(qobj, size, 10, 10);
		glTranslatef(-pos[i][0], -pos[i][1], -pos[i][2]);
	}	
}

template<typename T>
void TriDynamicMesh<T>::glDraw_StretchSpring(const RGBA& color, GLfloat width) const {
	Spring<T>::glDraw_Springs(stretchSpring, pos, color, width);
}

template<typename T>
void TriDynamicMesh<T>::glDraw_ShearSpring(const RGBA& color, GLfloat width) const {
	Spring<T>::glDraw_Springs(shearSpring, pos, color, width);
}

template<typename T>
void TriDynamicMesh<T>::glDraw_BendingSpring(const RGBA& color, GLfloat width) const {
	Spring<T>::glDraw_Springs(bendingSpring, pos, color, width);
}


};		// springspace

#endif