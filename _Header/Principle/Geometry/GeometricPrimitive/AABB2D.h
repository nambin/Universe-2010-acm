#ifndef _GEOMETRYSPACE_AABB2D_H_
#define _GEOMETRYSPACE_AABB2D_H_

#include "../../Arithmetic/Vector/Vector2.h"

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

namespace geometryspace {

using arithmetic::Vector2;

template<typename T>
class AABB2D {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructors
	inline AABB2D();
	template<typename R> explicit AABB2D(const AABB2D<R>&);
	inline AABB2D(const Vector2<T>& pt);
	inline AABB2D(const Vector2<T>& v1, const Vector2<T>& v2);

	//////////////////////////////////////////////////////////////////////////
	// Setter
	inline void init();
	inline void set(const Vector2<T>& pt);
	inline void set(const Vector2<T>& v1, const Vector2<T>& v2);

	//////////////////////////////////////////////////////////////////////////
	// Getter
	Vector2<T>& min()				{ return pp[0]; }
	Vector2<T>& max()				{ return pp[1]; }
	const Vector2<T>& min() const	{ return pp[0]; }
	const Vector2<T>& max() const	{ return pp[1]; }

	//////////////////////////////////////////////////////////////////////////
	// Manipulation
	inline void expand(T delta);
	inline void expand(const Vector2<T> &delta);

	inline void surround(const Vector2<T> &pt);
	inline void surround(const AABB2D &b);

	//////////////////////////////////////////////////////////////////////////
	// Utility
	inline const T area() const;
	inline const Vector2<T> center() const;
	inline const Vector2<T> extent() const;
	inline void boundingCircle(Vector2<T> &center, T &radius) const;	
	const int maximumExtent() const;

	inline void constrainedBdry(Vector2<T> &pt) const;

	//////////////////////////////////////////////////////////////////////////
	// Checker
	const bool isOverlaps(const AABB2D &b) const;
	const bool isInside(const Vector2<T>& pt) const;

	//////////////////////////////////////////////////////////////////////////
	// Friends Function
	template<typename T> friend const AABB2D surround(const AABB2D& b, const Vector2<T>& pt);
	template<typename T> friend const AABB2D surround(const AABB2D& b1, const AABB2D& b2);

public :
	Vector2<T> pp[2];
};

//////////////////////////////////////////////////////////////////////////
// Constructors
template<typename T> inline AABB2D<T>::AABB2D()												{ init();															}
template<typename T> template<typename R> AABB2D<T>::AABB2D(const AABB2D<R>& aabb)			{ pp[0] = Vector2<T>(aabb.pp[0]); pp[1] = Vector2<T>(aabb.pp[1]);	}
template<typename T> inline AABB2D<T>::AABB2D(const Vector2<T>& pt)							{ set(pt);															}
template<typename T> inline AABB2D<T>::AABB2D(const Vector2<T>& v1, const Vector2<T>& v2)	{ set(v1, v2);														}

template<typename T>
inline void AABB2D<T>::init() {
	pp[0] = Vector2<T>(static_cast<T>(basis::INFINITY));
	pp[1] = Vector2<T>(-static_cast<T>(basis::INFINITY));
}

template<typename T> inline void AABB2D<T>::set(const Vector2<T>& pt)						{ pp[0] = pp[1] = pt;							  }
template<typename T> inline void AABB2D<T>::set(const Vector2<T>& v1, const Vector2<T>& v2)	{ pp[0] = minVec(v1, v2); pp[1] = maxVec(v1, v2); }

template<typename T>
inline void AABB2D<T>::expand(T delta) {
	pp[0] -= Vector2<T>(delta);
	pp[1] += Vector2<T>(delta);
}

template<typename T>
inline void AABB2D<T>::expand(const Vector2<T> &delta)	{
	pp[0] -= delta;
	pp[1] += delta;
}

template<typename T> inline void AABB2D<T>::surround(const Vector2<T> &pt)	{ pp[0] = minVec(pp[0], pt);		pp[1] = maxVec(pp[1], pt);		}
template<typename T> inline void AABB2D<T>::surround(const AABB2D &b)		{ pp[0] = minVec(pp[0], b.min());	pp[1] = maxVec(pp[1], b.max());	}

template<typename T>
inline const T AABB2D<T>::area() const {
	Vector2<T> d = pp[1] - pp[0];
	return d.x() * d.y();
}

template<typename T>
inline const Vector2<T> AABB2D<T>::center() const {
	return static_cast<T>(0.5f) * (pp[0] + pp[1]);
}

template<typename T>
inline const Vector2<T> AABB2D<T>::extent() const {
	return pp[1] - pp[0];
}

template<typename T>
inline void AABB2D<T>::boundingCircle(Vector2<T> &center, T &radius) const {
	center = static_cast<T>(0.5f) * (pp[0] + pp[1]);
	radius = length(center, pp[0]);
}

template<typename T>
inline void AABB2D<T>::constrainedBdry(Vector2<T> &pt) const {
	for(Vector2<T>::index_type axis=0;axis<2;++axis) 
		pt[axis] = clamp(pt[axis], pp[0][axis], pp[1][axis]);
}

template<typename T>
const int AABB2D<T>::maximumExtent() const {
	Vector2<T> diag = pp[1] - pp[0];
	if(diag.x() > diag.y())	return 0;	
	else					return 1;
}

template<typename T>
const bool AABB2D<T>::isOverlaps(const AABB2D &b) const {
	bool x = (pp[1].x() >= b.pp[0].x()) && (pp[0].x() <= b.pp[1].x());
	bool y = (pp[1].y() >= b.pp[0].y()) && (pp[0].y() <= b.pp[1].y());	
	return (x && y);
}

template<typename T>
const bool AABB2D<T>::isInside(const Vector2<T>& pt)	const {
	return ((pp[0].x() <= pt.x() && pt.x() <= pp[1].x()) && 
			(pp[0].y() <= pt.y() && pt.y() <= pp[1].y())) 		
}

template<typename T>
inline const AABB2D<T> surround(const AABB2D<T> &b, const Vector2<T> &pt)	{ 
	return AABB2D<T>( minVec(b.min(), pt), maxVec(b.max(), pt));				
}

template<typename T>
inline const AABB2D<T> surround(const AABB2D<T>& b1, const AABB2D<T>& b2)	{ 
	return AABB2D<T>( minVec(b1.min(), b2.min()), maxVec(b1.max(), b2.max()));	
}

};		// geometryspace

#endif