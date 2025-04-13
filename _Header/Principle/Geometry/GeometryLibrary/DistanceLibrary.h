#ifndef _GEOMETRYSPACE_DISTANCE_LIBRARY_H_
#define _GEOMETRYSPACE_DISTANCE_LIBRARY_H_

#include "../../Principle.h"
#include "../../../Basis/Basis/BasicLibrary.h"
#include "../../Arithmetic/Vector/Vector2.h"
#include "../../Arithmetic/Vector/Vector3.h"

namespace geometryspace {

using arithmetic::Vector2;
using arithmetic::Vector3;

class PRINCIPLEDLL DistanceLibrary {
public : 
	template<typename T>
	static const T SqDistPointAABB(const Vector3<T> &p, const Vector3<T> &_min, const Vector3<T> &_max);

	template<typename T>
	static const T SqDistPointLine(const Vector2<T>& p, const Vector2<T>&, const Vector2<T>&);

enum Closest_Point_Type { CLOSEST_PT_SURFACE = -1, CLOSEST_PT_NODE0 = 0, CLOSEST_PT_NODE1 = 1, CLOSEST_PT_NODE2 = 2, CLOSEST_PT_EDGE0 =	3, CLOSEST_PT_EDGE1 = 4, CLOSEST_PT_EDGE2 =	5 };

	template<typename T>
	static const Vector3<T> ClosestPtPointTriangle(Vector3<T>& ret, Closest_Point_Type&, const Vector3<T> &p, const Vector3<T> &a, const Vector3<T> &b, const Vector3<T> &c);

	template<typename T>
	static const Vector3<T> ClosestPtPointTriangle(T ret[3], Closest_Point_Type&, const T p[3], const T a[3], const T b[3], const T c[3]);
};

template<typename T>
const T DistanceLibrary::SqDistPointAABB(const Vector3<T> &p, const Vector3<T> &_min, const Vector3<T> &_max) {
	T sqDist = T(0.0f);
	for(int i=0;i<3;++i) {
		T v = p[i];
		if(v < _min[i]) sqDist += (_min[i] - v) * (_min[i] - v);
		if(v > _max[i]) sqDist += (v - _max[i]) * (v - _max[i]);
	}
	return sqDist;
}

template<typename T>
const T DistanceLibrary::SqDistPointLine(const Vector2<T>& p, const Vector2<T>& v0, const Vector2<T>& v1) {
	static const T eps = 1.0e-8;

	const T Ai = v0[1] - v1[1];
	const T Bi = v1[0] - v0[0];
	const T Ci = - Ai * v0[0] - Bi * v0[1];

	const T eqn = Ai * p[0] + Bi * p[1] + Ci;

	if((SQ(v0[0]-v1[0]) + SQ(v0[1]-v1[1])) < eps)	{
//		std::cout << "DistanceLibrary::SqDistPointLine _ Segment is too short !!!!!!!!!!!!" << std::endl;
		return 1.0+6;
	}

	if(dot(p-v0, v1-v0) < 0)	return sqrlength(p, v0);
	if(dot(p-v1, v0-v1) < 0)	return sqrlength(p, v1);
	else						return SQ(eqn) / (SQ(Ai) + SQ(Bi));
}

template<typename T>
const Vector3<T> DistanceLibrary::ClosestPtPointTriangle(Vector3<T>& ret, Closest_Point_Type& pt_type, const Vector3<T> &p, const Vector3<T> &a, const Vector3<T> &b, const Vector3<T> &c ) {	
	return ClosestPtPointTriangle(ret.data, pt_type, p.data, a.data, b.data, c.data);	
}

template<typename T>
const Vector3<T> DistanceLibrary::ClosestPtPointTriangle( T result[3], Closest_Point_Type& closestPtType, const T p[3], const T a[3], const T b[3], const T c[3] )	{

	Vector3<T> bary;

	T ab[3];	SUB3(ab, b, a);
	T ac[3];	SUB3(ac, c, a);
	T ap[3];	SUB3(ap, p, a);

	T d1 = DOT3(ab, ap);
	T d2 = DOT3(ac, ap);
	if(d1 <= 0.0f && d2 <= 0.0f)	{	CPY3(result, a);		closestPtType = CLOSEST_PT_NODE0;		return Vector3<T>(1,0,0);	}

	T bp[3];	SUB3(bp, p, b);
	T d3 = DOT3(ab, bp);
	T d4 = DOT3(ac, bp);
	if(d3 >= 0.0f && d4 <= d3)		{	CPY3(result, b);		closestPtType = CLOSEST_PT_NODE1;		return Vector3<T>(0,1,0);	}

	T vc = d1 * d4 - d3 * d2;
	if(vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f)	{
		T v = d1 / (d1 - d3);
		MUL3(result, v, ab);		ADD3(result, result, a);
		closestPtType = CLOSEST_PT_EDGE2;	

		return Vector3<T>(1-v,v,0);
	}

	T cp[3];	SUB3(cp, p, c);
	T d5 = DOT3(ab, cp);
	T d6 = DOT3(ac, cp);
	if(d6 >= 0.0f && d5 <= d6)		{	CPY3(result, c);		closestPtType = CLOSEST_PT_NODE2;		return Vector3<T>(0,0,1);	}

	T vb = d5 * d2 - d1 * d6;
	if(vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f)	{
		T w = d2 / (d2 - d6);
		MUL3(result, w, ac);		ADD3(result, result, a);
		closestPtType = CLOSEST_PT_EDGE1;	

		return Vector3<T>(1-w,0,w);
	}

	T va = d3 * d6 - d5 * d4;
	if(va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f)	{
		T w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		T temp[3];	SUB3(temp, c, b);
		MUL3(result, w, temp);	ADD3(result, result, b);
		closestPtType = CLOSEST_PT_EDGE0;	

		return Vector3<T>(0,1-w,w);		
	}

	T denom = 1.0f / (va + vb + vc);
	T v = vb * denom;
	T w = vc * denom;
	MUL3(ab, v, ab);		MUL3(ac, w, ac);
	ADD3(result, a, ab);	ADD3(result, result, ac);
	closestPtType = CLOSEST_PT_SURFACE;	

	return Vector3<T>(1-v-w,v,w);
}

};		// namespace geometryspace

#endif