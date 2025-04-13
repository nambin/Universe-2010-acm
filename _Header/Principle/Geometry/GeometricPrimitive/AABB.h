#ifndef _GEOMETRYSPACE_AABB_H_
#define _GEOMETRYSPACE_AABB_H_

#include "../../Arithmetic/Vector/Vector3.h"

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

namespace geometryspace {

using arithmetic::Vector3;

template<typename T>
class AABB {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructors
	inline AABB();
	template<typename R> explicit AABB(const AABB<R>&);
	inline AABB(const Vector3<T>& pt);
	inline AABB(const Vector3<T>& v1, const Vector3<T>& v2);

	//////////////////////////////////////////////////////////////////////////
	// Setter
	inline void init();
	inline void set(const Vector3<T>& pt);
	inline void set(const Vector3<T>& v1, const Vector3<T>& v2);

	//////////////////////////////////////////////////////////////////////////
	// Getter
	Vector3<T>& min()				{ return pp[0]; }
	Vector3<T>& max()				{ return pp[1]; }
	const Vector3<T>& min() const	{ return pp[0]; }
	const Vector3<T>& max() const	{ return pp[1]; }

	//////////////////////////////////////////////////////////////////////////
	// Manipulation
	inline void expand(T delta);
	inline void expand(const Vector3<T> &delta);

	inline void surround(const Vector3<T> &pt);
	inline void surround(const AABB &b);

	//////////////////////////////////////////////////////////////////////////
	// Utility
	inline const T volume() const;
	inline const Vector3<T> center() const;
	inline const Vector3<T> extent() const;
	inline void boundingSphere(Vector3<T> &center, T &radius) const;	
	const int maximumExtent() const;

	inline void constrainedBdry(Vector3<T> &pt) const;

	//////////////////////////////////////////////////////////////////////////
	// Checker
	const bool isOverlaps(const AABB &b) const;
	const bool isInside(const Vector3<T>& pt) const;

	//////////////////////////////////////////////////////////////////////////
	// Friends Function
	template<typename T> friend const AABB surround(const AABB& b, const Vector3<T>& pt);
	template<typename T> friend const AABB surround(const AABB& b1, const AABB& b2);

public :
	Vector3<T> pp[2];
};

//////////////////////////////////////////////////////////////////////////
// Constructors
template<typename T> inline AABB<T>::AABB()												{ init();															}
template<typename T> template<typename R> AABB<T>::AABB(const AABB<R>& aabb)			{ pp[0] = Vector3<T>(aabb.pp[0]); pp[1] = Vector3<T>(aabb.pp[1]);	}
template<typename T> inline AABB<T>::AABB(const Vector3<T>& pt)							{ set(pt);															}
template<typename T> inline AABB<T>::AABB(const Vector3<T>& v1, const Vector3<T>& v2)	{ set(v1, v2);														}

template<typename T>
inline void AABB<T>::init() {
	pp[0] = Vector3<T>(static_cast<T>(basis::INFINITY));
	pp[1] = Vector3<T>(-static_cast<T>(basis::INFINITY));
}

template<typename T> inline void AABB<T>::set(const Vector3<T>& pt)							{ pp[0] = pp[1] = pt;							  }
template<typename T> inline void AABB<T>::set(const Vector3<T>& v1, const Vector3<T>& v2)	{ pp[0] = minVec(v1, v2); pp[1] = maxVec(v1, v2); }

template<typename T>
inline void AABB<T>::expand(T delta) {
	pp[0] -= Vector3<T>(delta);
	pp[1] += Vector3<T>(delta);
}

template<typename T>
inline void AABB<T>::expand(const Vector3<T> &delta)	{
	pp[0] -= delta;
	pp[1] += delta;
}

template<typename T> inline void AABB<T>::surround(const Vector3<T> &pt) { pp[0] = minVec(pp[0], pt);		pp[1] = maxVec(pp[1], pt);		}
template<typename T> inline void AABB<T>::surround(const AABB &b)		 { pp[0] = minVec(pp[0], b.min());	pp[1] = maxVec(pp[1], b.max());	}

template<typename T>
inline const T AABB<T>::volume() const {
	Vector3<T> d = pp[1] - pp[0];
	return d.x() * d.y() * d.z();
}

template<typename T>
inline const Vector3<T> AABB<T>::center() const {
	return static_cast<T>(0.5f) * (pp[0] + pp[1]);
}

template<typename T>
inline const Vector3<T> AABB<T>::extent() const {
	return pp[1] - pp[0];
}

template<typename T>
inline void AABB<T>::boundingSphere(Vector3<T> &center, T &radius) const {
	center = static_cast<T>(0.5f) * (pp[0] + pp[1]);
	radius = length(center, pp[0]);
}

template<typename T>
inline void AABB<T>::constrainedBdry(Vector3<T> &pt) const {
	for(Vector3<T>::index_type axis=0;axis<3;++axis) 
		pt[axis] = clamp(pt[axis], pp[0][axis], pp[1][axis]);
}

template<typename T>
const int AABB<T>::maximumExtent() const {
	Vector3<T> diag = pp[1] - pp[0];
	if(diag.x() > diag.y() && diag.x() > diag.z())	return 0;
	else if(diag.y() > diag.z())					return 1;
	else											return 2;
}

template<typename T>
const bool AABB<T>::isOverlaps(const AABB &b) const {
	bool x = (pp[1].x() >= b.pp[0].x()) && (pp[0].x() <= b.pp[1].x());
	bool y = (pp[1].y() >= b.pp[0].y()) && (pp[0].y() <= b.pp[1].y());
	bool z = (pp[1].z() >= b.pp[0].z()) && (pp[0].z() <= b.pp[1].z());
	return (x && y && z);
}

template<typename T>
const bool AABB<T>::isInside(const Vector3<T>& pt)	const {
	return ((pp[0].x() <= pt.x() && pt.x() <= pp[1].x()) && 
			(pp[0].y() <= pt.y() && pt.y() <= pp[1].y()) && 
			(pp[0].z() <= pt.z() && pt.z() <= pp[1].z()) );
}

template<typename T>
inline const AABB<T> surround(const AABB<T> &b, const Vector3<T> &pt)	{ 
	return AABB<T>( minVec(b.min(), pt), maxVec(b.max(), pt));				
}

template<typename T>
inline const AABB<T> surround(const AABB<T>& b1, const AABB<T>& b2)	{ 
	return AABB<T>( minVec(b1.min(), b2.min()), maxVec(b1.max(), b2.max()));	
}

};		// geometryspace

#endif