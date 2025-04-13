#pragma once

#include "../../Arithmetic/Vector/Vector3.h"
#include "GeometricObject.h"

namespace geometryspace {

using arithmetic::Vector3;

template<typename T>
class Cube : public GeometricObject<Vector3<T>> {
protected : 
	template<typename T> class CubeValue;

public : 
	Cube() : GeometricObject(new CubeValue<T>()) {}
	Cube(const Vector3<T>& _min, const Vector3<T>& _max) : GeometricObject(new CubeValue<T>(_min,_max)) {}

	static Material defaultMaterial();

protected :
	template<typename T>
	class CubeValue : public GeometricObjectValue<Vector3<T>> {
	public : 
		//////////////////////////////////////////////////////////////////////////
		// Constructors
		CubeValue() : GeometricObjectValue(Cube<T>::defaultMaterial()) {}
		CubeValue(const Vector3<T>& _min, const Vector3<T>& _max) : GeometricObjectValue(Cube<T>::defaultMaterial()), field(_min, _max) {}

		CubeValue * clone() const { return new CubeValue(*this); }

		//////////////////////////////////////////////////////////////////////////
		// Distance Field Library
		const bool			isInside(const Vector3<T>& p) const		{ return (getPhi(p) <= 0); }

		const T				getPhi(const Vector3<T>& p) const		{ return field(p); }
		const Vector3<T>	getNormal(const Vector3<T>& p) const	{ throwError("Not Implemented _ CubeValue::getNormal()"); return Vector3<T>(1,0,0); }

		//////////////////////////////////////////////////////////////////////////
		// glDraw
		void				glDraw() const;			

	public : 
		CubeDistanceField3D<T> field;
	};
};

template<typename T>
Material Cube<T>::defaultMaterial() {
	Material mat;
	mat.setAmbient(RGBA::red());
	mat.setDiffuse(RGBA::red());
	mat.setSpecular(RGBA::red());
	mat.setEmissive(RGBA::black());
	mat.setShininess(0);
	return mat;
}

template<typename R> template<typename T>
void Cube<R>::CubeValue<T>::glDraw() const {
	glEnable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_COLOR_MATERIAL);
	glShadeModel(GL_SMOOTH);

	mat.glBindMaterial();

	const Vector3<T> center = (field.get_min() + field.get_max()) * 0.5f;
	const Vector3<T> ext = field.get_max() - field.get_min();

	glTranslatef(center[0], center[1], center[2]);
	glScalef(ext[0], ext[1], ext[2]);
	glutSolidCube(1);
	glScalef(1.0f/ext[0], 1.0f/ext[1], 1.0f/ext[2]);
	glTranslatef(-center[0], -center[1], -center[2]);
}

};			// springspace

