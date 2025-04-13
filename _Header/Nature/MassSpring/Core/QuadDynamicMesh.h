#ifndef _MASS_SPRING_SPACE_QUAD_DYNAMICMESH_H_
#define _MASS_SPRING_SPACE_QUAD_DYNAMICMESH_H_

#include <vector>
#include "../../../Principle/Arithmetic/Transform/Transform.h"
#include "../Core/ClothLibrary.h"
#include "../Core/SpringLibrary.h"
#include "../Core/Constraints.h"
#include "../Core/MassSpringSystemMatrix.h"
#include "../Core/DynamicMesh.h"

namespace springspace {

using std::vector;
using basis::RGBA;
using basis::Array;
using arithmetic::Vector3;
using virtualspace::Material;

template<typename T>
class QuadDynamicMesh : public DynamicMesh {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	QuadDynamicMesh(index_type nX=0, index_type nY=0, T Lx=0, T Ly=0);
	~QuadDynamicMesh();

	void postProcessing() {}

	void init(index_type nX, index_type nY, T Lx, T Ly);
	void initpos();

	index_type dimX() const { return nX;	}
	index_type dimY() const { return nY;	}
	index_type size() const { return nX*nY; }

	float width() const	 { return (float)Lx; }
	float length() const { return (float)Ly; }

	index_type idx(index_type i, index_type j) const { return i * nY + j; }

	void setMass(T m)			{ mass = m / size();	}
	void setParticleMass(T m)	{ mass = m;				}

	//////////////////////////////////////////////////////////////////////////
	// Spring 
	void make_StretchSpring(SpringFunctor<T> *);
	void make_ShearSpring(SpringFunctor<T> *);
	void make_BendingSpring(SpringFunctor<T> *);

	void allocate_SystemMatrix(MassSpringSystemMatrix<T>&) const;

	//////////////////////////////////////////////////////////////////////////
	// Constraint
	void fixVertex(index_type i, index_type j);

	void setOnePointGrabConstraint();
	void setTwoPointGrabConstraint();

	//////////////////////////////////////////////////////////////////////////
	// Utility
	void transformMesh(const arithmetic::Transform<T> &trans);
	const Vector3<T> centerOfMass() const;

	const T maxStrain() const;

	void getVertexNormals(Array<Vector3<T>>&) const;

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	void glDraw_FlatSurface(const Material& mf, const Material& mb) const;
	void glDraw_SmoothSurface(const Material& mf, const Material& mb) const;

	void glDraw_Particle(const RGBA& c0, const RGBA& cConstraint, GLfloat size) const;
	void glDraw_StretchSpring(const RGBA&, GLfloat width) const;
	void glDraw_ShearSpring(const RGBA&, GLfloat width) const;
	void glDraw_BendingSpring(const RGBA&, GLfloat width) const;

public : 
	index_type					nX, nY;	
	T							dx, dy, dDiag;
	T							Lx, Ly;

	T							mass;

	Vec3Array<T>				pos, pos0;
	Vec3Array<T>				vel, vel0;

	SpringFunctor<T>*			stretchFunc;
	SpringFunctor<T>*			shearFunc; 
	SpringFunctor<T>*			bendingFunc;
	Array<Spring<T>>			stretchSpring;		
	Array<Spring<T>>			shearSpring;	
	Array<Spring<T>>			bendingSpring;	

	Constraints<T>				constraint;

	Array<vector<Spring<T>*>>	s2v;
};

template<typename T>
QuadDynamicMesh<T>::QuadDynamicMesh(index_type _nX, index_type _nY, T _Lx, T _Ly) 
: stretchFunc(NULL), shearFunc(NULL), bendingFunc(NULL)
{
	init(_nX, _nY, _Lx, _Ly);	
}

template<typename T>
QuadDynamicMesh<T>::~QuadDynamicMesh() {
	if(stretchFunc)	delete stretchFunc;	stretchFunc = NULL;
	if(shearFunc)	delete shearFunc;	shearFunc = NULL;
	if(bendingFunc)	delete bendingFunc;	bendingFunc = NULL;
}

template<typename T>
void QuadDynamicMesh<T>::init(index_type _nX, index_type _nY, T _Lx, T _Ly) {
	nX = _nX; nY = _nY; Lx = _Lx; Ly = _Ly; 
	
	if(nX <= 1 || nY <= 1)	{ dx = dy = 0;							}
	else					{ dx = Lx / (nX-1); dy = Ly / (nY-1);	}

	dDiag = sqrt(dx*dx + dy*dy);

	mass = 1.0f;

	pos.resize(nX*nY);	pos0.resize(nX*nY);
	vel.resize(nX*nY);	vel0.resize(nX*nY);
	initpos();

	if(stretchFunc)	delete stretchFunc;	stretchFunc = NULL;
	if(shearFunc)	delete shearFunc;	shearFunc = NULL;
	if(bendingFunc)	delete bendingFunc;	bendingFunc = NULL;

	stretchSpring.clear();
	shearSpring.clear();
	bendingSpring.clear();

	constraint.init(nX*nY);

	s2v.resize(nX*nY);
}

template<typename T>
void QuadDynamicMesh<T>::initpos() {
	for(QuadDynamicMesh<T>::index_type i=0;i<nX;++i)
		for(QuadDynamicMesh<T>::index_type j=0;j<nY;++j)
			pos(idx(i,j)) = pos0(idx(i,j)) = Vector3<T>((i*dx), (j*dy), 0);

	constraint.releaseConstraints();
}

template<typename T>
void QuadDynamicMesh<T>::make_StretchSpring(SpringFunctor<T> *func) {
	if(!stretchSpring.empty()) {
		std::cout << "Stretch Springs are already made" << std::endl;
		return;
	}

	stretchFunc = func;

	std::vector<Spring<T>> springs;
	for(QuadDynamicMesh<T>::index_type i=0;i<nX;++i) {
		for(QuadDynamicMesh<T>::index_type j=0;j<nY;++j) {			
			if(i != nX-1) springs.push_back(Spring<T>(idx(i,j), idx(i+1,j), dx));
			if(j != nY-1) springs.push_back(Spring<T>(idx(i,j), idx(i,j+1), dy));
		}
	}

	stretchSpring.resize(springs.size());
	stretchSpring.set(springs);	

	for(size_t i=0;i<stretchSpring.size();++i) {
		s2v[stretchSpring[i].idx0].push_back(&stretchSpring[i]);
		s2v[stretchSpring[i].idx1].push_back(&stretchSpring[i]);
	}
}

template<typename T>
void QuadDynamicMesh<T>::make_ShearSpring(SpringFunctor<T> *func) {
	if(!shearSpring.empty()) {
		std::cout << "Shear Springs are already made" << std::endl;
		return;
	}

	shearFunc = func;

	std::vector<Spring<T>> springs;
	for(QuadDynamicMesh<T>::index_type i=0;i<nX-1;++i) {
		for(QuadDynamicMesh<T>::index_type j=0;j<nY-1;++j) {			
			springs.push_back(Spring<T>(idx(i,j), idx(i+1,j+1), dDiag));
			springs.push_back(Spring<T>(idx(i+1,j), idx(i,j+1), dDiag));
		}
	}

	shearSpring.resize(springs.size());
	shearSpring.set(springs);	

	for(size_t i=0;i<shearSpring.size();++i) {
		s2v[shearSpring[i].idx0].push_back(&shearSpring[i]);
		s2v[shearSpring[i].idx1].push_back(&shearSpring[i]);
	}
}

template<typename T>
void QuadDynamicMesh<T>::make_BendingSpring(SpringFunctor<T> *func) {
}

template<typename T>
void QuadDynamicMesh<T>::allocate_SystemMatrix(MassSpringSystemMatrix<T>& mat) const {
	mat.clear();

	mat.n = nX * nY;
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
void QuadDynamicMesh<T>::fixVertex(index_type i, index_type j) {
	QuadDynamicMesh<T>::index_type _idx = idx(i,j);
	constraint.addConstraint(_idx, Vector3<T>(1,0,0));
	constraint.addConstraint(_idx, Vector3<T>(0,1,0));
	constraint.addConstraint(_idx, Vector3<T>(0,0,1));
	constraint.setVelocity(_idx, Vector3<T>(0));
	constraint.setHardConstraint(_idx);
}

template<typename T>
void QuadDynamicMesh<T>::setOnePointGrabConstraint() {
	fixVertex(0,0);
}

template<typename T>
void QuadDynamicMesh<T>::setTwoPointGrabConstraint() {
	fixVertex(0,0);
	fixVertex(nX-1,0);
}

template<typename T>
void QuadDynamicMesh<T>::transformMesh(const arithmetic::Transform<T> &trans) {
	for(size_t i=0;i<pos.size();++i)
		pos(i) = pos0(i) = trans.transformPoint(pos(i));
}

template<typename T>
const Vector3<T> QuadDynamicMesh<T>::centerOfMass() const {
	Vector3<T> com(0);
	for(size_t i=0;i<pos.size();++i)
		com += pos(i);
	com /= pos.size();
	return com;
}

template<typename T>
const T QuadDynamicMesh<T>::maxStrain() const {
	T _max = 0, strain;
	
	strain = SpringLibrary::maxStrain(stretchSpring, pos);	_max = (_max > strain) ? _max : strain;
	strain = SpringLibrary::maxStrain(shearSpring, pos);	_max = (_max > strain) ? _max : strain;
	strain = SpringLibrary::maxStrain(bendingSpring, pos);	_max = (_max > strain) ? _max : strain;

	return _max;
}

template<typename T>
void QuadDynamicMesh<T>::getVertexNormals(Array<Vector3<T>> &normals) const {
	assert(normals.size() == pos.size());

	const Vector3<T> *V0, *V1, *V2;
	Vector3<T> normal;	real area = 0;
	for(QuadDynamicMesh<T>::index_type i=0;i<nX;++i) {
		for(QuadDynamicMesh<T>::index_type j=0;j<nY;++j) {

			normals[idx(i,j)] = Vector3<T>(0);
			area = 0;

			if(i != nX-1 && j != nY-1) {
				V0 = &pos(idx(i,j));	V1 = &pos(idx(i+1,j));	V2 = &pos(idx(i,j+1));
				normal = cross(*V1-*V0, *V2-*V0);
				area += mag(normal);
				normals[idx(i,j)] += normal;
			}
			if(i != 0 && j != nY-1) {
				V0 = &pos(idx(i,j));	V1 = &pos(idx(i,j+1));	V2 = &pos(idx(i-1,j));
				normal = cross(*V1-*V0, *V2-*V0);
				area += mag(normal);
				normals[idx(i,j)] += normal;
			}
			if(i != 0 && j != 0) {
				V0 = &pos(idx(i,j));	V1 = &pos(idx(i-1,j));	V2 = &pos(idx(i,j-1));
				normal = cross(*V1-*V0, *V2-*V0);
				area += mag(normal);
				normals[idx(i,j)] += normal;
			}
			if(i != nX-1 && j != 0) {
				V0 = &pos(idx(i,j));	V1 = &pos(idx(i,j-1));	V2 = &pos(idx(i+1,j));
				normal = cross(*V1-*V0, *V2-*V0);
				area += mag(normal);
				normals[idx(i,j)] += normal;
			}

			normals[idx(i,j)] /= area;

		}
	}
}

template<typename T>
void QuadDynamicMesh<T>::glDraw_FlatSurface(const Material& mf, const Material& mb) const {	
	glEnable(GL_LIGHTING);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
	
	glDisable(GL_CULL_FACE);
	glDisable(GL_COLOR_MATERIAL);	

	mf.glBind_FrontMaterial();
	mb.glBind_BackMaterial();

	glBegin( GL_TRIANGLES );

	const Vector3<T> *V0, *V1, *V2;
	Vector3<T> normal;
	for(QuadDynamicMesh<T>::index_type i=0;i<nX-1;++i)	{
		for(QuadDynamicMesh<T>::index_type j=0;j<nY-1;++j)	{				
			V0 = &pos(idx(i,j));	V1 = &pos(idx(i+1,j));	V2 = &pos(idx(i,j+1));
			normal = unitVector(cross(*V1-*V0, *V2-*V0));

			glNormal3f(normal[0], normal[1], normal[2]);

			glVertex3f((*V0)[0], (*V0)[1], (*V0)[2]);
			glVertex3f((*V1)[0], (*V1)[1], (*V1)[2]);
			glVertex3f((*V2)[0], (*V2)[1], (*V2)[2]);

			V0 = &pos(idx(i+1,j));	V1 = &pos(idx(i+1,j+1));	V2 = &pos(idx(i,j+1));
			normal = unitVector(cross(*V1-*V0, *V2-*V0));

			glNormal3f(normal[0], normal[1], normal[2]);

			glVertex3f((*V0)[0], (*V0)[1], (*V0)[2]);
			glVertex3f((*V1)[0], (*V1)[1], (*V1)[2]);
			glVertex3f((*V2)[0], (*V2)[1], (*V2)[2]);
		}
	}

	glEnd();

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 0);
}

template<typename T>
void QuadDynamicMesh<T>::glDraw_SmoothSurface(const Material& mf, const Material& mb) const {
	Array<Vector3<T>> normals(pos.size());
	getVertexNormals(normals);

	glEnable(GL_LIGHTING);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);

	glDisable(GL_CULL_FACE);
	glDisable(GL_COLOR_MATERIAL);	

	mf.glBind_FrontMaterial();
	mb.glBind_BackMaterial();

	glBegin( GL_TRIANGLES );

	const Vector3<T> *V0, *V1, *V2;
	for(QuadDynamicMesh<T>::index_type i=0;i<nX-1;++i)	{
		for(QuadDynamicMesh<T>::index_type j=0;j<nY-1;++j)	{				
			V0 = &pos(idx(i,j));	V1 = &pos(idx(i+1,j));	V2 = &pos(idx(i,j+1));

			glNormal3f(normals[idx(i,j)][0], normals[idx(i,j)][1], normals[idx(i,j)][2]);
			glVertex3f((*V0)[0], (*V0)[1], (*V0)[2]);
			glNormal3f(normals[idx(i+1,j)][0], normals[idx(i+1,j)][1], normals[idx(i+1,j)][2]);
			glVertex3f((*V1)[0], (*V1)[1], (*V1)[2]);
			glNormal3f(normals[idx(i,j+1)][0], normals[idx(i,j+1)][1], normals[idx(i,j+1)][2]);
			glVertex3f((*V2)[0], (*V2)[1], (*V2)[2]);

			V0 = &pos(idx(i+1,j));	V1 = &pos(idx(i+1,j+1));	V2 = &pos(idx(i,j+1));

			glNormal3f(normals[idx(i+1,j)][0], normals[idx(i+1,j)][1], normals[idx(i+1,j)][2]);
			glVertex3f((*V0)[0], (*V0)[1], (*V0)[2]);
			glNormal3f(normals[idx(i+1,j+1)][0], normals[idx(i+1,j+1)][1], normals[idx(i+1,j+1)][2]);
			glVertex3f((*V1)[0], (*V1)[1], (*V1)[2]);
			glNormal3f(normals[idx(i,j+1)][0], normals[idx(i,j+1)][1], normals[idx(i,j+1)][2]);
			glVertex3f((*V2)[0], (*V2)[1], (*V2)[2]);
		}
	}

	glEnd();

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 0);
}

template<typename T>
void QuadDynamicMesh<T>::glDraw_Particle(const RGBA& c0, const RGBA& c1, GLfloat size) const {
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
void QuadDynamicMesh<T>::glDraw_StretchSpring(const RGBA& color, GLfloat width) const {
	Spring<T>::glDraw_Springs(stretchSpring, pos, color, width);
}

template<typename T>
void QuadDynamicMesh<T>::glDraw_ShearSpring(const RGBA& color, GLfloat width) const {
	Spring<T>::glDraw_Springs(shearSpring, pos, color, width);
}

template<typename T>
void QuadDynamicMesh<T>::glDraw_BendingSpring(const RGBA& color, GLfloat width) const {
	Spring<T>::glDraw_Springs(bendingSpring, pos, color, width);
}

};		// springspace

#endif