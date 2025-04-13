#ifndef _ARITHMETIC_TRANSFORMS_H_
#define _ARITHMETIC_TRANSFORMS_H_

#include <cmath>
#include "../../../Basis/MemoryManagement/RCPtr.h"
#include "../Matrix/Matrix4.h"

namespace arithmetic {

using basis::RCObject;
using basis::RCPtr;
using arithmetic::Vector3;
using arithmetic::Vector4;
using arithmetic::Matrix4;

template<typename T> class Transform {
public :
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	inline Transform();
	inline Transform(const Matrix4<T> &mat);
	inline Transform(const T mat[4][4]);
	inline Transform(const Matrix4<T>& mat, const Matrix4<T>& mInv);

	//////////////////////////////////////////////////////////////////////////
	// Setter methods
	inline void set(const Transform& t);
	inline void set(const Matrix4<T> &mat);
	inline void set(const T mat[4][4]);
	inline void set(const Matrix4<T>& mat, const Matrix4<T>& mInv);

	inline void identity();

	//////////////////////////////////////////////////////////////////////////
	// Getter methods
	inline const Matrix4<T>* const getMatrixPtr() const		{ return &value->m;		}
	inline const Matrix4<T>* const getMatrixInvPtr() const	{ return &value->mInv;	}

	inline const Matrix4<T>& getMatrixRef() const		{ return value->m;		}
	inline const Matrix4<T>& getMatrixInvRef() const	{ return value->mInv;	}

	//////////////////////////////////////////////////////////////////////////
	// Special Transform Setting
	inline const Transform getInverse() const; 
	inline const Transform &translate(const Vector3<T> &v);
	inline const Transform &scale(T x, T y, T z);
	inline const Transform &rotateX(T angle);
	inline const Transform &rotateY(T angle);
	inline const Transform &rotateZ(T angle);
	inline const Transform &rotate(T angle, const Vector3<T> &axis);
	const Transform &lookAt(const Vector3<T> &pos, const Vector3<T> &look, const Vector3<T> &up);

	//////////////////////////////////////////////////////////////////////////
	// Transform Operation
	const Vector3<T> transformPoint(const Vector3<T> &pt) const;
	const Vector3<T> transformVector(const Vector3<T> &v) const;
	const Vector3<T> transformNormal(const Vector3<T> &normal) const;
	
	const Transform operator*(const Transform &t) const;

private : 
	template<typename T>
	struct TransformValue : public RCObject {
		TransformValue()							{ m.identity(); mInv.identity();			}
		TransformValue(const TransformValue &rhs)	{ m = rhs.m;	mInv = rhs.mInv;			}
		TransformValue(const Matrix4<T> &mat)		{ m = mat;		mInv = m.getInverse();		}
		TransformValue(const T mat[4][4])			{ m = mat;		mInv = m.getInverse();		}
		TransformValue(const Matrix4<T> &mat, const Matrix4<T> &matInv) { m = mat;	mInv = matInv;	}

		TransformValue *clone() const { return new TransformValue(*this); }

		~TransformValue() {}

		Matrix4<T> m, mInv;
	};

	//////////////////////////////////////////////////////////////////////////
	// Reference Counted Pointer
	RCPtr<TransformValue<T>> value;

public : 
	const RCPtr<TransformValue<T>>& getTransformValue() const { return value; }
};

template<typename T>
inline Transform<T>::Transform() : value(new Transform<T>::TransformValue<T>()) {}

template<typename T>
inline Transform<T>::Transform(const Matrix4<T> &mat) : value(new Transform<T>::TransformValue<T>(mat)) {}

template<typename T>
inline Transform<T>::Transform(const T mat[4][4]) : value(new Transform<T>::TransformValue<T>(mat)) {}

template<typename T>
inline Transform<T>::Transform(const Matrix4<T>& mat, const Matrix4<T>& matInv) : value(new Transform<T>::TransformValue<T>(mat, matInv)) {}

template<typename T>
inline void Transform<T>::set(const Transform<T>& t) { 
	(*this) = t; 
}

template<typename T>
inline void Transform<T>::set(const Matrix4<T> &mat) { 
	value = new Transform::TransformValue<T>(mat);
}

template<typename T>
inline void Transform<T>::set(const T mat[4][4]) { 
	value = new Transform::TransformValue<T>(mat);
}

template<typename T>
inline void Transform<T>::set(const Matrix4<T>& mat, const Matrix4<T>& mInv) { 
	value = new Transform::TransformValue<T>(mat, mInv);
}

template<typename T>
inline void Transform<T>::identity() { 
	value = new Transform::TransformValue<T>();
}

template<typename T>
inline const Transform<T> Transform<T>::getInverse() const {
	Transform t(value->mInv, value->m);
	return t;
}

template<typename T>
inline const Transform<T> &Transform<T>::translate(const Vector3<T> &v) {
	Transform t;
	t.value->m.setTranslation(v);
	t.value->mInv.setTranslation(-v);
	return ((*this)=t);
}

template<typename T>
inline const Transform<T> &Transform<T>::scale(T x, T y, T z) {
	Transform t;
	t.value->m.setScale(Vector3<T>(x,y,z));
	t.value->mInv.setScale(Vector3<T>(1.0f/x, 1.0f/y, 1.0f/z));
	return ((*this)=t);
}

template<typename T>
inline const Transform<T> &Transform<T>::rotateX(T angle)	{
	Transform t;
	t.value->m.setRotationX(angle);
	t.value->mInv = t.value->m.getTranspose();
	return ((*this)=t);
}

template<typename T>
inline const Transform<T> &Transform<T>::rotateY(T angle)	{
	Transform t;
	t.value->m.setRotationY(angle);
	t.value->mInv = t.value->m.getTranspose();
	return ((*this)=t);
}

template<typename T>
inline const Transform<T> &Transform<T>::rotateZ(T angle)	{
	Transform t;
	t.value->m.setRotationZ(angle);
	t.value->mInv = t.value->m.getTranspose();
	return ((*this)=t);
}

template<typename T>
inline const Transform<T> &Transform<T>::rotate(T angle, const Vector3<T> &axis)	{
	Transform t;
	t.value->m.setRotationAxis(angle, axis);
	t.value->mInv = t.value->m.getTranspose();
	return ((*this)=t);
}

template<typename T>
const Transform<T> &Transform<T>::lookAt(const Vector3<T> &pos, const Vector3<T> &look, const Vector3<T> &up)	{

	Transform t;
	Matrix4<T> &m = t.value->m;
	Matrix4<T> &mInv = t.value->mInv;

	m.setView(pos, look, up);

	mInv.set(0,3,pos.x());	mInv.set(1,3,pos.y());	mInv.set(2,3,pos.z());	mInv.set(3,3,1.0f);

	Vector3<T> zAxis = unitVector(pos - look);
	Vector3<T> xAxis = unitVector(cross(up, zAxis));
	Vector3<T> yAxis = cross(zAxis, xAxis);

	mInv.set(0,0,xAxis.x());	mInv.set(1,0,xAxis.y());	mInv.set(2,0,xAxis.z());	mInv.set(3,0,0.0f);
	mInv.set(0,1,yAxis.x());	mInv.set(1,1,yAxis.y());	mInv.set(2,1,yAxis.z());	mInv.set(3,1,0.0f);
	mInv.set(0,2,zAxis.x());	mInv.set(1,2,zAxis.y());	mInv.set(2,2,zAxis.z());	mInv.set(3,2,0.0f);

	return ((*this)=t);
}

template<typename T>
const Vector3<T> Transform<T>::transformPoint(const Vector3<T> &pt) const	{
	Vector4<T> _pt(pt, 1.0f);
	Vector4<T> _ptTrans = value->m * _pt;
	if(_ptTrans[3] == 1.0f)	return _ptTrans.xyz();
	else					return (_ptTrans.xyz()/_ptTrans[3]);
}

template<typename T>
const Vector3<T> Transform<T>::transformVector(const Vector3<T> &v) const	{
	Vector4<T> _v(v, 0.0f);
	Vector4<T> _vTrans = value->m * _v;
	return _vTrans.xyz();
}

template<typename T>
const Vector3<T> Transform<T>::transformNormal(const Vector3<T> &normal) const	{

	Matrix4<T> &mInv = value->mInv;

	T x = normal.x(), y = normal.y(), z = normal.z();
	return Vector3<T>(mInv.get(0,0) * x + mInv.get(1,0) * y + mInv.get(2,0) * z,
					  mInv.get(0,1) * x + mInv.get(1,1) * y + mInv.get(2,1) * z,
					  mInv.get(0,2) * x + mInv.get(1,2) * y + mInv.get(2,2) * z);
}

template<typename T>
inline const Transform<T> Transform<T>::operator*(const Transform &t) const {
	Matrix4<T> m1 = value->m * t.value->m;
	Matrix4<T> m2 = t.value->mInv * value->mInv;

	return Transform(m1, m2);
}

template<typename T>
inline const Transform<T> translate(const Vector3<T> &v)	{
	Transform<T> trans;
	trans.translate(v);
	return trans;	
}

template<typename T>
inline const Transform<T> scale(T x, T y, T z)	{
	Transform<T> trans;
	trans.scale(x,y,z);
	return trans;	
}

template<typename T>
inline const Transform<T> rotateX(T angle)	{
	Transform<T> trans;
	trans.rotateX(angle);
	return trans;	
}

template<typename T>
inline const Transform<T> rotateY(T angle)	{
	Transform<T> trans;
	trans.rotateY(angle);
	return trans;	
}

template<typename T>
inline const Transform<T> rotateZ(T angle)	{
	Transform<T> trans;
	trans.rotateZ(angle);
	return trans;	
}

template<typename T>
inline const Transform<T> rotate(T angle, const Vector3<T> &axis)	{
	Transform<T> trans;
	trans.rotate(angle, axis);
	return trans;	
}

template<typename T>
inline const Transform<T> lookAt(const Vector3<T> &pos, const Vector3<T> &look, const Vector3<T> &up)	{
	Transform<T> trans;
	trans.lookAt(pos, look, up);
	return trans;	
}

template<typename T>
const Transform<T> orthographic(T zNear, T zFar) {
	const T _one = static_cast<T>(1.0f);
	const T _zero = static_cast<T>(0.0f);
	return scale(_one, _one, 2.0f / (zNear-zFar)) * translate(Vector3<T>(_zero, _zero, (zNear + zFar) * 0.5f));
}

template<typename T>
const Transform<T> perspective(T fov, T n, T f) {
	T inv_denom = 1.0f / (n-f);
	Matrix4<T> persp(1.0f,0.0f,0.0f,0.0f,0.0f,1.0f,0.0f,0.0f,0.0f,0.0f,(n+f)*inv_denom,-1.0f,0.0f,0.0f,2.0f*f*n*inv_denom,0);

	T invTanAng = 1.0f / tan(Rad(fov) * 0.5f);

	return scale(invTanAng, invTanAng, real(1.0f)) * Transform<T>(persp);	
}

};		// geometryspace

#endif