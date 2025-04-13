#ifndef _GEOMETRYSPACE_GEOMETRY_SPHERE_H_
#define _GEOMETRYSPACE_GEOMETRY_SPHERE_H_
	
#include "../../Arithmetic/Vector/Vector3.h"
#include "GeometricObject.h"

namespace geometryspace {

using arithmetic::Vector3;

template<typename T>
class Sphere : public GeometricObject<Vector3<T>> {
protected : 
	template<typename T> class SphereValue;

public : 
	Sphere() : GeometricObject(new SphereValue<T>()) {}
	Sphere(const Vector3<T>& c, T r) : GeometricObject(new SphereValue<T>(c,r)) {}

	static Material defaultMaterial();

protected :
	template<typename T>
	class SphereValue : public GeometricObjectValue<Vector3<T>> {
	public : 
		//////////////////////////////////////////////////////////////////////////
		// Constructors
		SphereValue() : GeometricObjectValue(Sphere<T>::defaultMaterial()), radius(1) {}
		SphereValue(const Vector3<T>& c, T r) : GeometricObjectValue(Sphere<T>::defaultMaterial()), center(c), radius(r) {}

		SphereValue * clone() const { return new SphereValue(*this); }

		//////////////////////////////////////////////////////////////////////////
		// Distance Field Library
		const bool			isInside(const Vector3<T>& p) const		{ return (getPhi(p) <= 0); }

		const T				getPhi(const Vector3<T>& p) const		{ return arithmetic::length(p, center) - radius; }
		const Vector3<T>	getNormal(const Vector3<T>& p) const	{ return unitVector(p - center); }

		//////////////////////////////////////////////////////////////////////////
		// glDraw
		void				glDraw() const;			

	public : 
		Vector3<T>	center;
		T			radius;
	};
};

template<typename T>
Material Sphere<T>::defaultMaterial() {
	Material mat;
	mat.setAmbient(RGBA::red());
	mat.setDiffuse(RGBA::red());
	mat.setSpecular(RGBA::red());
	mat.setEmissive(RGBA::black());
	mat.setShininess(0);
	return mat;
}

template<typename R> template<typename T>
void Sphere<R>::SphereValue<T>::glDraw() const {
	glEnable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_COLOR_MATERIAL);
	glShadeModel(GL_SMOOTH);
	
	mat.glBindMaterial();

	GLUquadricObj *qobj = gluNewQuadric();

	gluQuadricDrawStyle(qobj, (GLenum) GLU_FILL);
	gluQuadricNormals(qobj, (GLenum) GLU_SMOOTH);

	glTranslatef(center[0], center[1], center[2]);		
	gluSphere(qobj, radius, 30, 30);
	glTranslatef(-center[0], -center[1], -center[2]);
}

};			// springspace

#endif