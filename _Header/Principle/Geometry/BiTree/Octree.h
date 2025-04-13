#ifndef _GEOMETRY_OCTREE_H_
#define _GEOMETRY_OCTREE_H_

#include "OctreeCell.h"

namespace geometryspace {

#define DEFAULT_MAX_DEPTH 10

//////////////////////////////////////////////////////////////////////////
// Octree Base Class
template<typename _REAL, typename _CELL>
class _Octree {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	_Octree(unsigned int ml, Vector3<_REAL> lc, Vector3<_REAL> uc);
	virtual ~_Octree();

	void init();
	void init(unsigned int ml, Vector3<_REAL> lc, Vector3<_REAL> uc);
	void release();

	//////////////////////////////////////////////////////////////////////////
	// Utility
	void constrainedBdry(Vector3<_REAL> &pos) const;

	//////////////////////////////////////////////////////////////////////////
	// Get Location Code _ pos should lie between lc & uc
	inline Vector3<_REAL> getPos(const unsigned int loc[3]) const;
	inline void getLocCode(unsigned int loc[3], const Vector3<_REAL> &pos) const;

	inline Vector3<_REAL> getCellSize(unsigned int level) const;

	inline Vector3<_REAL> getCenterPos(const _CELL *cell) const;

	//////////////////////////////////////////////////////////////////////////
	// Locate Function _ pos, p0, p1 should lie in the domain
	inline const _CELL * getLocateCell(const Vector3<_REAL> &pos) const;
	const _CELL * getLocateRegion(const Vector3<_REAL> &p0, const Vector3<_REAL> &p1) const;

	//////////////////////////////////////////////////////////////////////////
	// Make Tree
	template<typename SPLITPROC> void makeTree(const SPLITPROC &sproc);
	template<typename SPLITPROC> void splitCell(_CELL *, const SPLITPROC &sproc);

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	void glDrawOctree(const RGBA &color, GLfloat width = 1.0f) const;
	void glDrawOnlyMinimumLeaf(const RGBA &color, GLfloat width = 1.0f) const;
	virtual void glDraw(const RGBA& color) const {}

protected : 
	//////////////////////////////////////////////////////////////////////////
	// Basic Split Procedure
	class _SplitProc {};

	template<typename _OCTREE>
	class SplitProc : public _SplitProc {
	public : 
		SplitProc(_OCTREE *_oct) : oct(_oct) {}

		_CELL * createRootCell(unsigned int rootLevel) const;		

		inline const bool isSplit(const _CELL *cell) const;

		void splitCell(_CELL *cell) const;
		void postProcess() const;		

	public : 
		_OCTREE *oct;
	};

public : 
	virtual _SplitProc* getSplitProc() const = 0;

public : 	
	unsigned int			maxLevel, rootLevel;
	_REAL					maxVal;

	Vector3<_REAL>			lc, uc, range;

	size_t					numCells;
	_CELL*					root;
	std::vector<_CELL*>		cells;
};

//////////////////////////////////////////////////////////////////////////
// Constructor/Destructor
template<typename _REAL, typename _CELL>
_Octree<_REAL, _CELL>::_Octree(unsigned int ml, Vector3<_REAL> lc, Vector3<_REAL> uc) : root(NULL) {	
	init(ml, lc, uc);	
}

template<typename _REAL, typename _CELL>
_Octree<_REAL, _CELL>::~_Octree() {
	release();
}

template<typename _REAL, typename _CELL>
void _Octree<_REAL, _CELL>::init() {
	init(maxLevel, lc, uc);
}

template<typename _REAL, typename _CELL>
void _Octree<_REAL, _CELL>::init(unsigned int _ml, Vector3<_REAL> _lc, Vector3<_REAL> _uc) {
	release();
	maxLevel = _ml;
	rootLevel = _ml - 1;
	maxVal = pow(2.0f, (int)rootLevel);
	lc = _lc;
	uc = _uc; 
	range = uc-lc;
	root = NULL;
}

template<typename _REAL, typename _CELL>
void _Octree<_REAL, _CELL>::release() {
	cells.clear();

	if(root)
		geometryspace::release(root);	
	root = NULL;
}

//////////////////////////////////////////////////////////////////////////
// Utility
template<typename _REAL, typename _CELL>
void _Octree<_REAL, _CELL>::constrainedBdry(Vector3<_REAL> &pos) const {
	if ( pos[0] < lc[0] ) pos[0] = lc[0]; else if ( pos[0] > uc[0] ) pos[0] = uc[0];
	if ( pos[1] < lc[1] ) pos[1] = lc[1]; else if ( pos[1] > uc[1] ) pos[1] = uc[1];
	if ( pos[2] < lc[2] ) pos[2] = lc[2]; else if ( pos[2] > uc[2] ) pos[2] = uc[2];
}

//////////////////////////////////////////////////////////////////////////
// Get Location Code _ pos should lie between lc & uc
template<typename _REAL, typename _CELL>
inline Vector3<_REAL> _Octree<_REAL, _CELL>::getPos(const unsigned int loc[3]) const {	
	Vector3<_REAL> p(loc[0] / maxVal, loc[1] / maxVal, loc[2] / maxVal);
	return (lc + range * p);							// does need Range Checker ?
}

template<typename _REAL, typename _CELL>
inline void _Octree<_REAL, _CELL>::getLocCode(unsigned int loc[3], const Vector3<_REAL> &pos) const {
	Vector3<_REAL> p = (pos - lc) / range;

	for(size_t i=0;i<3;++i) {
		loc[i] = (unsigned int)(p[i] * maxVal);
		if(loc[i] & (1 << rootLevel))
			--loc[i];
	}	
}

template<typename _REAL, typename _CELL>
inline Vector3<_REAL> _Octree<_REAL, _CELL>::getCellSize(unsigned int level) const {
	unsigned int divide = (1 << (rootLevel - level));
	return range / static_cast<_REAL>(divide);
}

template<typename _REAL, typename _CELL>
inline Vector3<_REAL> _Octree<_REAL, _CELL>::getCenterPos(const _CELL *cell) const {
	return getPos(cell->locCode) + getCellSize(cell->level) * static_cast<_REAL>(0.5f);
}

//////////////////////////////////////////////////////////////////////////
// Locate Function _ pos, p0, p1 should lie in the domain
template<typename _REAL, typename _CELL>
inline const _CELL * _Octree<_REAL, _CELL>::getLocateCell(const Vector3<_REAL> &pos) const {
	unsigned int loc[3];
	getLocCode(loc, pos);
	return OctreeLibrary::getLocateCell(root, loc);
}

template<typename _REAL, typename _CELL>
const _CELL * _Octree<_REAL, _CELL>::getLocateRegion(const Vector3<_REAL> &p0, const Vector3<_REAL> &p1) const {
	unsigned int loc0[3], loc1[3];
	getLocCode(loc0, p0);
	getLocCode(loc1, p1);

	unsigned int diff[3] = { loc0[0] ^ loc1[0], loc0[1] ^ loc1[1], loc0[2] ^ loc1[2] };

	int level = rootLevel;
	int minLevel = rootLevel;
	while ((level != -1) && !(diff[0] & (1 << level))) level--;
	while ((minLevel > level) && !(diff[1] & (1 << minLevel))) minLevel--;

	return OctreeLibrary::getLocateCell(root, loc0, minLevel + 1);
}	

//////////////////////////////////////////////////////////////////////////
// Make Tree
template<typename _REAL, typename _CELL> template<typename SPLITPROC>
void _Octree<_REAL, _CELL>::makeTree(const SPLITPROC &sproc) {
	release();

	root = sproc.createRootCell(rootLevel);	
	cells.push_back(root);

	unsigned int idx = 0;
	while(idx < cells.size()) {
		splitCell(cells[idx], sproc);
		++idx;
	}

	sproc.postProcess();
}

template<typename _REAL, typename _CELL> template<typename SPLITPROC>
void _Octree<_REAL, _CELL>::splitCell(_CELL *cell, const SPLITPROC &sproc) {
	if(!sproc.isSplit(cell))
		return;

	sproc.splitCell(cell);
}

//////////////////////////////////////////////////////////////////////////
// glDraw
template<typename _REAL, typename _CELL> 
void _Octree<_REAL, _CELL>::glDrawOctree(const RGBA &color, GLfloat width) const{
	for(unsigned int i=0;i<cells.size();++i) {
		if(OctreeLibrary::isLeaf(cells[i])) {
			OctreeLibrary::glDrawCell(cells[i], this, color, width);
		}
	}
}

template<typename _REAL, typename _CELL> 
void _Octree<_REAL, _CELL>::glDrawOnlyMinimumLeaf(const RGBA &color, GLfloat width) const{
	for(unsigned int i=0;i<cells.size();++i) {
		if(cells[i]->level == 0) {
			OctreeLibrary::glDrawCell(cells[i], this, color, width);
		}
	}
}

//////////////////////////////////////////////////////////////////////////
// Basic Split Procedure
template<typename _REAL, typename _CELL> template<typename _OCTREE>
_CELL * _Octree<_REAL, _CELL>::SplitProc<_OCTREE>::createRootCell(unsigned int rootLevel) const {
	unsigned int rootCode[3] = {0,0,0};
	return new _CELL(rootCode, rootLevel, NULL);	
}

template<typename _REAL, typename _CELL> template<typename _OCTREE>
inline const bool _Octree<_REAL, _CELL>::SplitProc<_OCTREE>::isSplit(const _CELL *cell) const {
	return !(cell->level == 0);
	//if(cell->level == 0)
	//	return false;
	//return (rand() / (double)RAND_MAX) < 0.4f;
}

template<typename _REAL, typename _CELL> template<typename _OCTREE>
void _Octree<_REAL, _CELL>::SplitProc<_OCTREE>::splitCell(_CELL *cell) const { 
	cell->children = new _CELL[8];
	for(unsigned int i=0;i<8;++i) {
		cell->children[i].parent = cell;
		cell->children[i].children = NULL;
		cell->children[i].level = cell->level - 1;

		for(int axis=0;axis<3;++axis) 
			cell->children[i].locCode[axis] = (cell->locCode[axis] | (((i & (1 << axis)) >> axis) << cell->children[i].level));

		oct->cells.push_back(&cell->children[i]);
	}
}

template<typename _REAL, typename _CELL> template<typename _OCTREE>
void _Octree<_REAL, _CELL>::SplitProc<_OCTREE>::postProcess() const {
	oct->numCells = oct->cells.size();
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Base Octree /////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
template<typename _REAL>
class Octree : public _Octree<_REAL, OctCell> {
public : 
	Octree(unsigned int ml = 2, Vector3<_REAL> lc = Vector3<_REAL>(0), Vector3<_REAL> uc = Vector3<_REAL>(1))
	: _Octree(ml, lc, uc) {}

public : 
	typename _Octree<_REAL, OctCell>::SplitProc<Octree> * getSplitProc() const;	
};

template<typename _REAL>
typename _Octree<_REAL, OctCell>::SplitProc<Octree<_REAL>> * Octree<_REAL>::getSplitProc() const {
	return new _Octree<_REAL, OctCell>::SplitProc<Octree<_REAL>>(const_cast<Octree<_REAL>*>(this));
}

};

#endif