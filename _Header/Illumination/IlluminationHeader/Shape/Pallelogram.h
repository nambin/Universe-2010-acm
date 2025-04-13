#ifndef _ILLUMINATION_PARALLELOGRAM_H_
#define _ILLUMINATION_PARALLELOGRAM_H_

#include "Shape.h"

namespace illuminationspace {

class ILLUMINATIONDLL Parallelogram : public Shape {
public :
	static inline const real ParallelEpsilon()	{ return 1.0e-3f;	}
	static inline const real Epsilon()			{ return 1.0e-3f;	}	

public : 
	Parallelogram(const Transform &o2w, const Vector3<real>& _base, const Vector3<real>& _u, const Vector3<real>& _v);
	Parallelogram(const Transform &o2w, const Vector3<real>& _base, const Vector3<real>& _u, const Vector3<real>& _v, 
				  const Vector2<real> & _uv0, const Vector2<real>& _uv1, const Vector2<real>& _uv2);

protected :
	class ILLUMINATIONDLL ParallelogramValue : public Shape::ShapeValue {
	public : 
		ParallelogramValue(const Transform &o2w, const Vector3<real>& _base, const Vector3<real>& _u, const Vector3<real>& _v);
		ParallelogramValue(const Transform &o2w, const Vector3<real>& _base, const Vector3<real>& _u, const Vector3<real>& _v, 
						   const Vector2<real> & _uv0, const Vector2<real>& _uv1, const Vector2<real>& _uv2);

		ParallelogramValue *clone() const { return new ParallelogramValue(*this); }

		const BBox objectBound() const;
		const BBox worldBound() const;

		const real area() const { return _area; }	
		const real pdf() const	{ return _pdf;	}

		const Vector3<real> sample(real u1, real u2, Vector3<real> *N) const;
		const ONB &sample(real u1, real u2, Vector3<real> &pos) const;

		const bool hit(const Ray& r, real &tHit, DifferentialGeometry &geometry) const;
		const bool shadowHit(const Ray& r) const;

	public : 
		Vector3<real>	base;		//
		Vector3<real>	u;			// One Vertex and the offsets to its neighbors
		Vector3<real>	v;			//

		Vector3<real>	norm;		// normal of the parallelogram
		Vector3<real>	unorm;		// normalized u offset
		Vector3<real>	vnorm;		// normalized v offset

		ONB				uvw;

		Vector2<real>	uv0;		//
		Vector2<real>	uv1;		// texture coords
		Vector2<real>	uv2;		//

		real			_area;
		real			_pdf;		// Precomputed const pdf
	};
};

};		// illuminationspace

#endif
