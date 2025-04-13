#ifndef _GEOMETRYSPACE_GEOMETRY_POLE_H_
#define _GEOMETRYSPACE_GEOMETRY_POLE_H_

#include "../../Arithmetic/Vector/Vector3.h"
#include "GeometricObject.h"

namespace geometryspace {

using arithmetic::Vector3;

template<typename T>
class Cylinder : public GeometricObject<Vector3<T>> {
protected : 
	template<typename T> class CylinderValue;

public : 
	Cylinder() : GeometricObject(new CylinderValue<T>()) {}
	Cylinder(T radius, T length) : GeometricObject(new CylinderValue<T>(radius,length)) {}

	void set_rigid_transformation(const arithmetic::Transform<T>& transform) { static_cast<RCPtr<CylinderValue<T>>>(value)->set_rigid_transformation(transform); }

	static Material defaultMaterial();

protected :
	template<typename T>
	class CylinderValue : public GeometricObjectValue<Vector3<T>> {
	public : 
		//////////////////////////////////////////////////////////////////////////
		// Constructors
		CylinderValue() : GeometricObjectValue(Cylinder<T>::defaultMaterial()), field(1,1) {}
		CylinderValue(T _radius, T _length) : GeometricObjectValue(Cylinder<T>::defaultMaterial()), field(_radius, _length) {}

		CylinderValue * clone() const { return new CylinderValue(*this); }

		void set_rigid_transformation(const arithmetic::Transform<T>& transform) { field.set_rigid_transformation(transform); }

		//////////////////////////////////////////////////////////////////////////
		// Distance Field Library
		const bool			isInside(const Vector3<T>& p) const		{ return (getPhi(p) <= 0);	}

		const T				getPhi(const Vector3<T>& p) const		{ return field.getPhi(p);	}
		const Vector3<T>	getNormal(const Vector3<T>& p) const	{ throwError("Not Implemented _ CylinderValue::getNormal()"); return Vector3<T>(1,0,0);	}

		//////////////////////////////////////////////////////////////////////////
		// glDraw
		void				glDraw() const;			

	public : 
		CylinderDistanceField<T>	field;
	};
};

template<typename T>
Material Cylinder<T>::defaultMaterial() {
	Material mat;
	mat.setAmbient(RGBA::red());
	mat.setDiffuse(RGBA::red());
	mat.setSpecular(RGBA::red());
	mat.setEmissive(RGBA::black());
	mat.setShininess(0);
	return mat;
}

template<typename R> template<typename T>
void Cylinder<R>::CylinderValue<T>::glDraw() const {

	GLfloat matrix[16];
	for(int i=0;i<16;++i) matrix[i] = field.get_rigid_transform().getMatrixRef().getptr()[i];

	glEnable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_COLOR_MATERIAL);
	glShadeModel(GL_SMOOTH);

	mat.glBindMaterial();

	GLUquadricObj *qobj = gluNewQuadric();

	gluQuadricDrawStyle(qobj, (GLenum) GLU_FILL);
	gluQuadricNormals(qobj, (GLenum) GLU_SMOOTH);

	glMatrixMode(GL_MODELVIEW);

	glPushMatrix();

	glMultMatrixf(matrix);	

	glTranslatef(0, -field.get_length() * 0.5f, 0);
	glRotatef(-90, 1, 0, 0);
	
	glPushMatrix();	glRotatef(180, 1, 0, 0); gluDisk(qobj, 0, field.get_radius(), 32, 1); glPopMatrix();
	glPushMatrix();	glTranslatef(0, 0, field.get_length()); gluDisk(qobj, 0, field.get_radius(), 32, 1); glPopMatrix();

	gluCylinder(qobj, field.get_radius(), field.get_radius(), field.get_length(), 32, 1);	
	
	glPopMatrix();
}

};			// springspace

#endif