#ifndef _GEOMETRYSPACE_PLANE_H_
#define _GEOMETRYSPACE_PLANE_H_

#include "../../Arithmetic/Vector/Vector3.h"
#include "../GeometryLibrary/ONB.h"
#include "GeometricObject.h"

namespace geometryspace {

using arithmetic::Vector3;

template<typename T>
class Plane : public GeometricObject<Vector3<T>> {
protected : 
	template<typename T> class PlaneValue;

public : 
	Plane() : GeometricObject(new PlaneValue<T>()) {}
	Plane(const Vector3<T>& n, T d) : GeometricObject(new PlaneValue<T>(n,d)) {}

	static Material defaultFrontMaterial();
	static Material defaultBackMaterial();

protected :
	template<typename T>
	class PlaneValue : public GeometricObjectValue<Vector3<T>> {
	public : 
		//////////////////////////////////////////////////////////////////////////
		// Constructors
		PlaneValue() 
		: GeometricObjectValue(Plane<T>::defaultFrontMaterial()), backMaterial(Plane<T>::defaultBackMaterial()), normal(0,1,0), dist(0) {}

		PlaneValue(const Vector3<T>& n, T d) 
		: GeometricObjectValue(Plane<T>::defaultFrontMaterial()), backMaterial(Plane<T>::defaultBackMaterial()), normal(arithmetic::unitVector(n)), dist(d / arithmetic::mag(n)) {}

		PlaneValue *		clone() const { return new PlaneValue(*this); }

		//////////////////////////////////////////////////////////////////////////
		// Distance Field Library
		const bool			isInside(const Vector3<T>& p) const		{ return (getPhi(p) <= 0);					}

		const T				getPhi(const Vector3<T>& p) const		{ return arithmetic::dot(p, normal) - dist;	}
		const Vector3<T>	getNormal(const Vector3<T>& p) const	{ return normal;							}

		//////////////////////////////////////////////////////////////////////////
		// glDraw
		void				glDraw() const;	
	
	public : 
		Vector3<T>	normal;	
		T			dist;

		Material	backMaterial;
	};
};

template<typename T>
Material Plane<T>::defaultFrontMaterial() {
	Material mat;
	mat.setAmbient(RGBA::darkBrown());
	mat.setDiffuse(RGBA::darkBrown());
	mat.setSpecular(RGBA::darkBrown());
	mat.setEmissive(RGBA::black());
	mat.setShininess(0);
	return mat;
}

template<typename T>
Material Plane<T>::defaultBackMaterial() {
	Material mat;
	mat.setAmbient(RGBA::paleBrown());
	mat.setDiffuse(RGBA::paleBrown());
	mat.setSpecular(RGBA::paleBrown());
	mat.setEmissive(RGBA::black());
	mat.setShininess(0);
	return mat;
}

template<typename R> template<typename T>
void Plane<R>::PlaneValue<T>::glDraw() const {
	glEnable(GL_LIGHTING);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);

	glEnable(GL_DEPTH_TEST);

	glDisable(GL_CULL_FACE);
	glDisable(GL_COLOR_MATERIAL);
	glShadeModel(GL_SMOOTH);

	mat.glBind_FrontMaterial();
	backMaterial.glBind_BackMaterial();

	T c = 1.0e+3;

	ONB<T> uvw;
	uvw.initFromUnitVector(normal);

	Vector3<T> ori = normal * dist;
	Vector3<T> p0 = ori - c * uvw.U - c * uvw.V;
	Vector3<T> p1 = ori + c * uvw.U - c * uvw.V;
	Vector3<T> p2 = ori + c * uvw.U + c * uvw.V;
	Vector3<T> p3 = ori - c * uvw.U + c * uvw.V;

	glNormal3f(normal[0], normal[1], normal[2]);
	glBegin(GL_QUADS);
		glVertex3f(p0[0], p0[1], p0[2]);
		glVertex3f(p1[0], p1[1], p1[2]);
		glVertex3f(p2[0], p2[1], p2[2]);
		glVertex3f(p3[0], p3[1], p3[2]);
	glEnd();

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 0);
}

};			// springspace

#endif