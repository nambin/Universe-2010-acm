#ifndef _GEOMETRY_BOUNDING_VOLUME_HIERARCHY_H_
#define _GEOMETRY_BOUNDING_VOLUME_HIERARCHY_H_

#include <vector>
#include "../GeometricPrimitive/AABB.h"

namespace geometryspace {

using std::vector;

template<typename DATA, typename T>
class BVH {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructors
	BVH();	
	BVH(BVH* lt, BVH* rt, bool dataCopy = false);
	BVH(const AABB<T>&);
	virtual ~BVH();

	void	release();
	void	set(BVH* lt, BVH* rt, bool dataCopy = false);

	//////////////////////////////////////////////////////////////////////////
	// Getter
	AABB<T>&			getBBox()		{ return bbox; }
	const AABB<T>&		getBBox() const { return bbox; }

	BVH*				getLeftBVH()		{ return left; }
	const BVH*			getLeftBVH() const	{ return left; }
	BVH*				getRightBVH()		{ return right;}
	const BVH*			getRightBVH() const	{ return right;}

	vector<DATA>&		getData()		{ return data; }
	const vector<DATA>& getData() const { return data; }

	//////////////////////////////////////////////////////////////////////////
	// Setter
	void				pushData(const DATA&);
	void				pushData(const vector<DATA>&);

	//////////////////////////////////////////////////////////////////////////
	// Query
	bool				isRoot() const { (parent == NULL);					}
	bool				isLeaf() const { (left == NULL && right == NULL);	}

	//////////////////////////////////////////////////////////////////////////
	// Manipulate
	void				removeDuplicateData();

protected :
	BVH*			parent;
	BVH*			left;
	BVH*			right;
	
	AABB<T>			bbox;
	vector<DATA>	data;
};

template<typename DATA, typename T>
BVH<DATA,T>::BVH() : parent(NULL), left(NULL), right(NULL) {}

template<typename DATA, typename T>
BVH<DATA,T>::BVH(BVH* lt, BVH* rt, bool dataCopy) : parent(NULL), left(NULL), right(NULL) {
	set(lt, rt, dataCopy);
}

template<typename DATA, typename T>
BVH<DATA,T>::BVH(const AABB<T>& _bbox) : parent(NULL), left(NULL), right(NULL) {
	bbox = _bbox;
}

template<typename DATA, typename T>
BVH<DATA,T>::~BVH() {
	release();
}

template<typename DATA, typename T>
void BVH<DATA,T>::release() {
	parent = NULL;
	if(left) delete left;
	if(right) delete right;
	
	left = right = NULL;

	bbox.init();
	data.clear();
}

template<typename DATA, typename T>
void BVH<DATA,T>::set(BVH* lt, BVH* rt, bool dataCopy = false) {
	release();

	bbox = geometryspace::surround(lt->bbox, rt->bbox);
	left = lt;
	right = rt;

	if(left)	left->parent = this;
	if(right)	right->parent = this;

	if(dataCopy) {
		if(left)	pushData(left->data);
		if(right)	pushData(right->data);
	}	
}

template<typename DATA, typename T>
void BVH<DATA,T>::pushData(const DATA& da) {
	data.push_back(da);
}

template<typename DATA, typename T>
void BVH<DATA,T>::pushData(const vector<DATA>& da) {
	for(size_t i=0;i<da.size();++i)
		pushData(da[i]);
}

template<typename DATA, typename T>
void BVH<DATA,T>::removeDuplicateData() {
	data.clear();
	if(left)  left->removeDuplicateData();
	if(right) right->removeDuplicateData();
}

};		// geometryspace

#endif
