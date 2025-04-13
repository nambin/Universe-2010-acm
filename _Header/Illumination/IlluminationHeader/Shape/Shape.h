#ifndef _ILLUMINATION_SHAPE_H_
#define _ILLUMINATION_SHAPE_H_

#include <vector>
#include "../../../Basis/MemoryManagement/RCPtr.h"
#include "../../../Principle/Arithmetic/Vector/Vector2.h"
#include "../../../Principle/Arithmetic/Vector/Vector3.h"
#include "../Illumination/Random.h"
#include "../Core/BBox.h"
#include "../Core/ONB.h"
#include "../Core/Transform.h"

namespace illuminationspace {

using std::vector;
using arithmetic::Vector2;
using arithmetic::Vector3;

class Ray;

struct DifferentialGeometry	{
	Vector3<real>	hit_p;
	ONB				uvw;
	Vector3<real>	shadingNormal;

	Vector2<real>	uv;
};

class ILLUMINATIONDLL Shape {
protected : 
	class ShapeValue;

public :
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Shape() : value(NULL) {}	
	Shape(ShapeValue *ptr) : value(ptr) {}

	//////////////////////////////////////////////////////////////////////////
	// Bounding Box
	const BBox objectBound() const	{ return value->objectBound();	}
	const BBox worldBound() const	{ return value->worldBound();	}

	//////////////////////////////////////////////////////////////////////////
	// Manipulator	
	void refine(std::vector<Shape> &refined) const	{ value->refine(refined);		}

	void getShadingGeometry(const Transform &obj2world, const DifferentialGeometry &dg, DifferentialGeometry &dgShading) const
	{ value->getShadingGeometry(obj2world, dg, dgShading); }

	//////////////////////////////////////////////////////////////////////////
	// Property
	const bool canIntersect() const	{ return value->canIntersect();	}
	const real area() const			{ return value->area();	}
	const real pdf() const			{ return value->pdf();	}

	//////////////////////////////////////////////////////////////////////////
	// Sampling
	const Vector3<real> sample(real u1, real u2, Vector3<real> *N) const							{ return value->sample(u1, u2, N);		}
	const Vector3<real> sample(const Vector3<real> &pt, real u1, real u2, Vector3<real> *N) const	{ return value->sample(pt, u1, u2, N);	}

	const ONB& sample(real u1, real u2, Vector3<real> &pos) const	{ return value->sample(u1, u2, pos);	}

	//////////////////////////////////////////////////////////////////////////
	// Intersection Test
	const bool hit(const Ray& r, real& tHit, DifferentialGeometry& dg) const	{ return value->hit(r, tHit, dg);	}
	const bool shadowHit(const Ray& r) const									{ return value->shadowHit(r);		}

protected :
	class ILLUMINATIONDLL ShapeValue : public basis::RCObject {
	public : 
		ShapeValue() : ObjectToWorld(Transform()), WorldToObject(Transform()) {}
		ShapeValue(const Transform &o2w);
		virtual ~ShapeValue() {}

		virtual ShapeValue *clone() const = 0;		

		//////////////////////////////////////////////////////////////////////////
		// Bounding Box
		virtual const BBox objectBound() const = 0;
		virtual const BBox worldBound() const;

		//////////////////////////////////////////////////////////////////////////
		// Manipulator		
		virtual void refine(std::vector<Shape> &refined) const;
		virtual void getShadingGeometry(const Transform &obj2world, const DifferentialGeometry &dg, DifferentialGeometry &dgShading) const;

		//////////////////////////////////////////////////////////////////////////
		// Property
		virtual const bool canIntersect() const	{ return true;	}
		virtual const real area() const;
		virtual const real pdf() const;

		//////////////////////////////////////////////////////////////////////////
		// Sampling
		virtual const Vector3<real> sample(real u1, real u2, Vector3<real> *N) const;							// used for sampling Area Light
		virtual const Vector3<real> sample(const Vector3<real> &pt, real u1, real u2, Vector3<real> *N) const;	// Need to check _ When is used ?

		virtual const ONB& sample(real u1, real u2, Vector3<real> &pos) const;									// used for generating photon path

		//////////////////////////////////////////////////////////////////////////
		// Intersection Test
		virtual const bool hit(const Ray& r, real& tHit, DifferentialGeometry& geometry) const;
		virtual const bool shadowHit(const Ray& r) const;		

	public : 
		Transform ObjectToWorld, WorldToObject;
	};	

	//////////////////////////////////////////////////////////////////////////
	// Reference Counted Pointer
	basis::RCPtr<ShapeValue> value;

public : 
	const basis::RCPtr<ShapeValue> &getShapeValue() const { return value; }	
};

class ILLUMINATIONDLL ShapeSet : public Shape {
protected : 
	class ShapeSetValue;

public : 	
	ShapeSet(ShapeSetValue* ptr = NULL) : Shape(ptr) {}
	ShapeSet(const std::vector<Shape> &s, const Transform &o2w)	: Shape(new ShapeSetValue(s, o2w)) {}

protected : 
	class ILLUMINATIONDLL ShapeSetValue : public Shape::ShapeValue {
	public:	
		//////////////////////////////////////////////////////////////////////////
		// Constructor/Destructor
		ShapeSetValue();
		ShapeSetValue(const Transform &o2w);
		ShapeSetValue(const std::vector<Shape> &s, const Transform &o2w);

		ShapeSetValue * clone() const { return new ShapeSetValue(*this); }

	protected : 
		void calcAreas();

	public : 
		//////////////////////////////////////////////////////////////////////////
		// Bounding Box
		const BBox objectBound() const;

		//////////////////////////////////////////////////////////////////////////
		// Manipulator		
		void refine(std::vector<Shape> &refined) const;

		//////////////////////////////////////////////////////////////////////////
		// Property
		const bool canIntersect() const;
		const real area() const	{ return _area;	}

		//////////////////////////////////////////////////////////////////////////
		// Sampling
		const Vector3<real> sample(real u1, real u2, Vector3<real> *Ns) const;			// used for sampling Area Light

		//////////////////////////////////////////////////////////////////////////
		// Intersection Test
		const bool hit(const Ray &ray, real &t_hit, DifferentialGeometry &dg) const;
		const bool shadowHit(const Ray &ray) const;

	public :
		// ShapeSet Private Data
		real			_area;
		vector<real>	areaCDF;
		vector<Shape>	shapes;		
	};
};

};		// illuminationspace

#endif
