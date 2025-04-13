#ifndef _MASS_SPRING_SPACE_SPRING__H_
#define _MASS_SPRING_SPACE_SPRING__H_

#include <gl/GL.h>
#include "../../Nature.h"
#include "../../../Basis/Basis/Array.h"
#include "../../../Basis/Color/RGBA.h"
#include "../../../Principle/Arithmetic/Vector/Vector3.h"

namespace springspace {

using basis::RGBA;
using basis::Array;
using arithmetic::Vector3;

template<typename T>
class Spring {
public : 
	Spring() {}
	Spring(unsigned int idx0, unsigned int idx1, T restlen);

	static 
	void glDraw_Springs(const Array<Spring<T>>& springs, const Vec3Array<T>& pos, const RGBA &, GLfloat width);

public : 
	unsigned int	idx0, idx1;
	T				restlen;	
};

template<typename T>
inline Spring<T>::Spring(unsigned int _i0, unsigned int _i1, T _restlen) 
: idx0(_i0), idx1(_i1), restlen(_restlen) {}

template<typename T>
void Spring<T>::glDraw_Springs(const Array<Spring<T>>& springs, const Vec3Array<T>& pos, const RGBA& color, GLfloat width) {	
	glDisable(GL_LIGHTING);
		
	glColor4f( color[0], color[1], color[2], color[3] );
	glLineWidth(width);

	glBegin(GL_LINES);
	for(size_t i=0;i<springs.size();++i) {
		const Vector3<T>& v0 = pos[springs[i].idx0];
		const Vector3<T>& v1 = pos[springs[i].idx1];

		glVertex3f(v0[0], v0[1], v0[2]);
		glVertex3f(v1[0], v1[1], v1[2]);
	}
	glEnd();
}

};		// springspace

#endif