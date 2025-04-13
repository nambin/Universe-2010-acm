#ifndef _GEOMETRYSPACE_GEOMETRY_CIRCLE_H_
#define _GEOMETRYSPACE_GEOMETRY_CIRCLE_H_

#include "../../Arithmetic/Vector/Vector2.h"
#include "GeometricObject.h"

namespace geometryspace {

using arithmetic::Vector2;

template<typename T>
class Circle : public GeometricObject<Vector2<T>> {
protected : 
	template<typename T> class CircleValue;

public : 
	Circle() : GeometricObject(new CircleValue<T>()) {}
	Circle(const Vector2<T>& c, T r) : GeometricObject(new CircleValue<T>(c,r)) {}

	static Material defaultMaterial();

protected :
	template<typename T>
	class CircleValue : public GeometricObjectValue<Vector2<T>> {
	public : 
		//////////////////////////////////////////////////////////////////////////
		// Constructors
		CircleValue() : GeometricObjectValue(Circle<T>::defaultMaterial()), radius(1) {}
		CircleValue(const Vector2<T>& c, T r) : GeometricObjectValue(Circle<T>::defaultMaterial()), center(c), radius(r) {}

		CircleValue * clone() const { return new CircleValue(*this); }

		//////////////////////////////////////////////////////////////////////////
		// Distance Field Library
		const bool			isInside(const Vector2<T>& p) const		{ return (getPhi(p) <= 0); }

		const T				getPhi(const Vector2<T>& p) const		{ return arithmetic::length(p, center) - radius; }
		const Vector2<T>	getNormal(const Vector2<T>& p) const	{ return unitVector(p - center); }

		//////////////////////////////////////////////////////////////////////////
		// glDraw
		void				glDraw() const;			

	public : 
		Vector2<T>	center;
		T			radius;
	};
};

template<typename T>
Material Circle<T>::defaultMaterial() {
	Material mat;
	mat.setAmbient(RGBA::brown());
	mat.setDiffuse(RGBA::brown());
	mat.setSpecular(RGBA::brown());
	mat.setEmissive(RGBA::black());
	mat.setShininess(0);
	return mat;
}

template<typename R> template<typename T>
void Circle<R>::CircleValue<T>::glDraw() const {

	const unsigned int n_divide = 40;

	const RGBA c = mat.getDiffuse();
	glColor4f(c[0], c[1], c[2], c[3]);

	glBegin(GL_POLYGON);
	for(unsigned int i=0;i<n_divide;++i) {
		const T r0 = T(i) / T(n_divide);
		const Vector2<T> p = center + radius * Vector2<T>(cos(TWOPI * r0), sin(TWOPI * r0));

		glVertex3f(p[0], p[1], 0.001f);
	}
	glEnd();
}

};			// springspace

#endif