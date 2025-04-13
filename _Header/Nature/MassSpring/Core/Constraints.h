#ifndef _MASS_SPRING_SPACE_CONSTRAINTS_H_
#define _MASS_SPRING_SPACE_CONSTRAINTS_H_

#include <vector>
#include "../../../Basis/Basis/Array2.h"
#include "../Core/Vec3Array.h"

namespace springspace {

using std::vector;
using basis::Array2;

template<typename T>
class Constraints {
public : 
	Constraints(size_t size = 0);

	void init(size_t size);
	void clear();

	//////////////////////////////////////////////////////////////////////////
	// Manipulate Constraints
	void releaseConstraints();
	void releaseSoftConstraints();

	void releaseConstraint(size_t i);	

	void addConstraint(size_t i, Vector3<T> dir);
	void setVelocity(size_t i, const Vector3<T>& vel);
	void setHardConstraint(size_t i) { isHardConstraint[i] = true; }

	//////////////////////////////////////////////////////////////////////////
	// Use Constraints
	const short			numConstraintDirection(size_t i) const	{ return constrainDOF[i];			}
	const bool			isConstrained(size_t i) const			{ return (constrainDOF[i] != 0);	}

	const Vector3<T>&	getVelocity(size_t i) const;	
	const Vector3<T>	getVelocity(size_t i, const Vector3<T>& vel) const;	

	void				filterVelocities(Vec3Array<T>& vel) const;	
	void				filterVelocity(size_t i, Vector3<T>& vel) const;	

	void				update_Constrained_IdxList() const;

private : 
	Array<bool>				isHardConstraint;
	Array<short>			constrainDOF;
	Array2<Vector3<T>>		dir;
	Vec3Array<T>			vel;

	// Auxiliary Variable for Velocity Filtering
	mutable bool			isValid_IdxList;		// flag for checking validity of constrained_IdxList[]
	mutable vector<size_t>	constrained_IdxList;
};

template<typename T>
Constraints<T>::Constraints(size_t size) {
	init(size);
}

template<typename T>
void Constraints<T>::init(size_t size) {
	isHardConstraint.resize(size);
	constrainDOF.resize(size);
	dir.resize(size, 3);
	vel.resize(size);

	for(size_t i=0;i<size;++i) {
		isHardConstraint[i] = false;
		constrainDOF[i] = 0;
		vel[i] = Vector3<T>(0);
	}

	isValid_IdxList = true;
	constrained_IdxList.clear();
}

template<typename T>
void Constraints<T>::clear() {
	init(0);
}

//////////////////////////////////////////////////////////////////////////
// Manipulate Constraints
template<typename T>
void Constraints<T>::releaseConstraints() {
	for(size_t i=0;i<constrainDOF.size();++i)
		releaseConstraint(i);
}

template<typename T>
void Constraints<T>::releaseSoftConstraints() {
	for(size_t i=0;i<constrainDOF.size();++i) {
		if(!isHardConstraint[i]) {
			releaseConstraint(i);
		}
	}
}

template<typename T>
void Constraints<T>::releaseConstraint(size_t idx) {
	isHardConstraint[idx] = false;
	constrainDOF[idx] = 0;
	vel[idx] = Vector3<T>(0);

	// Update Auxiliary Variable
	isValid_IdxList = false;
}

template<typename T>
void Constraints<T>::addConstraint(size_t idx, Vector3<T> _dir) {
	if(constrainDOF[idx] > 2)
		return;

	for(short i=0;i<constrainDOF[idx];++i) 
		_dir -= (dot(_dir, dir(idx, i)) * dir(idx, i));

	if(sqrmag(_dir) < 1.0e-4)
		return;

	dir(idx, constrainDOF[idx]++) = unitVector(_dir);

	// Update Auxiliary Variable
	isValid_IdxList = false;
}

template<typename T>
void Constraints<T>::setVelocity(size_t idx, const Vector3<T>& _vel) {
	vel[idx] = Vector3<T>(0);
	for(short i=0;i<constrainDOF[idx];++i) 
		vel[idx] += (dot(_vel, dir(idx, i)) * dir(idx, i));
}

//////////////////////////////////////////////////////////////////////////
// Use Constraints
template<typename T>
const Vector3<T>& Constraints<T>::getVelocity(size_t idx) const {
	return vel[idx];
}

template<typename T>
const Vector3<T> Constraints<T>::getVelocity(size_t idx, const Vector3<T> &_vel) const {
	if(constrainDOF[idx] == 0)
		return _vel;

	Vector3<T> v(_vel);
	filterVelocity(idx, v);
	return v + vel[idx];
}

template<typename T>
void Constraints<T>::filterVelocities(Vec3Array<T>& _vel) const {
	if(!isValid_IdxList)
		update_Constrained_IdxList();

	for(size_t i=0;i<constrained_IdxList.size();++i)
		filterVelocity(constrained_IdxList[i], _vel[constrained_IdxList[i]]);
}

template<typename T>
void Constraints<T>::filterVelocity(size_t idx, Vector3<T> &_vel) const {
	for(short i=0;i<constrainDOF[idx];++i) 
		_vel -= (dot(_vel, dir(idx, i)) * dir(idx, i));
}

template<typename T>
void Constraints<T>::update_Constrained_IdxList() const {
	if(isValid_IdxList)
		return;

	constrained_IdxList.clear();
	for(size_t i=0;i<constrainDOF.size();++i) {
		if(isConstrained(i)) {
			constrained_IdxList.push_back(i);
		}
	}

	isValid_IdxList = true;
}

};		// springspace

#endif
