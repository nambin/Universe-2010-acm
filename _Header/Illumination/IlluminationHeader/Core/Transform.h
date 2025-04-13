#ifndef _ILLUMINATION_TRANSFORMS_H_
#define _ILLUMINATION_TRANSFORMS_H_

#include <cmath>
#include "../../../Principle/Arithmetic/Transform/Transform.h"
#include "Ray.h"
#include "BBox.h"
#include "ONB.h"
#include "../../IlluminationSpace.h"
#include "../Illumination/Constants.h"

namespace illuminationspace {

using arithmetic::Matrix4;

class ILLUMINATIONDLL Transform : public arithmetic::Transform<real> {
public :
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	inline Transform();
	inline Transform(const Transform &trans);
	inline Transform(const arithmetic::Transform<real> &trans);
	inline Transform(const Matrix4<real> &mat);
	inline Transform(const real mat[4][4]);
	inline Transform(const Matrix4<real>& mat, const Matrix4<real>& matInv);

	//////////////////////////////////////////////////////////////////////////
	// Transformer
	const Ray transformRay(const Ray& r) const;
	const BBox transformBBox(const BBox &bbox) const;
	const ONB transformONB(const ONB &uvw) const;

	void transformRay(const Ray& r, Ray* rTrans) const;

	inline Transform &operator=(const arithmetic::Transform<real> &t);

#ifdef _USE_XML_
	//////////////////////////////////////////////////////////////////////////
	// XML IO
	static const std::string XmlElementName() { return std::string("Transform"); }
	void importXmlNode(xmlNodePtr node);

private : 
	void getTransformFromXmlNode(xmlNodePtr node);

#endif
};

inline Transform::Transform() : arithmetic::Transform<real>() {}
inline Transform::Transform(const Transform &trans) : arithmetic::Transform<real>(trans) {}
inline Transform::Transform(const arithmetic::Transform<real> &trans) : arithmetic::Transform<real>(trans) {}
inline Transform::Transform(const Matrix4<real> &mat) : arithmetic::Transform<real>(mat) {}
inline Transform::Transform(const real mat[4][4]) : arithmetic::Transform<real>(mat) {}
inline Transform::Transform(const Matrix4<real>& mat, const Matrix4<real>& matInv) : arithmetic::Transform<real>(mat,matInv) {}

inline Transform &Transform::operator=(const arithmetic::Transform<real> &t) {
	arithmetic::Transform<real>::set(t);
	return (*this);
}

};		// illuminationspace

#endif
