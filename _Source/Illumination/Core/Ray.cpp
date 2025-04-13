#include <stdafx.h>
#include <Illumination/IlluminationHeader/Core/Ray.h>

using namespace illuminationspace;

void Ray::setDirection(const Vector3<real>& v)	{
	data[1] = v;
	data[2] = Vector3<real>(  (v.x() == 0) ? 1.0e+6 : 1.0f / v.x()
							, (v.y() == 0) ? 1.0e+6 : 1.0f / v.y()
							, (v.z() == 0) ? 1.0e+6 : 1.0f / v.z() );

	posneg[0] = (data[1].x() > 0 ? 0 : 1);
	posneg[1] = (data[1].y() > 0 ? 0 : 1);
	posneg[2] = (data[1].z() > 0 ? 0 : 1);
}