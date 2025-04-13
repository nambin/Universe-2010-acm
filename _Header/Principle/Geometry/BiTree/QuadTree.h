#ifndef _GEOMETRY_QUADTREE_H_
#define _GEOMETRY_QUADTREE_H_

#include "QuadTreeCell.h"

namespace geometryspace {

#define DEFAULT_MAX_DEPTH 10

//////////////////////////////////////////////////////////////////////////
// QuadTree Base Class
template<typename _REAL, typename _CELL>
class _QuadTree {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	_QuadTree(unsigned int ml, Vector2<_REAL> lc, Vector2<_REAL> uc);
	virtual ~_QuadTree();

	void init();
	void init(unsigned int ml, Vector2<_REAL> lc, Vector2<_REAL> uc);
	void release();

	//////////////////////////////////////////////////////////////////////////
	// Utility
	void constrainedBdry(Vector2<_REAL> &pos) const;

	//////////////////////////////////////////////////////////////////////////
	// Get Location Code _ pos should lie between lc & uc
	inline Vector2<_REAL> getPos(const unsigned int loc[2]) const;
	inline void getLocCode(unsigned int loc[2], const Vector2<_REAL> &pos) const;

	inline Vector2<_REAL> getCellSize(unsigned int level) const;

	inline Vector2<_REAL> getCenterPos(const _CELL *cell) const;
	
	//////////////////////////////////////////////////////////////////////////
	// Locate Function _ pos, p0, p1 should lie in the domain
	inline const _CELL * getLocateCell(const Vector2<_REAL> &pos) const;
	const _CELL * getLocateRegion(const Vector2<_REAL> &p0, const Vector2<_REAL> &p1) const;

	//////////////////////////////////////////////////////////////////////////
	// Make Tree
	template<typename SPLITPROC> void makeTree(const SPLITPROC &sproc);
	template<typename SPLITPROC> void splitCell(_CELL *, const SPLITPROC &sproc);

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	void glDrawQuadTree(const RGBA &color, GLfloat width = 1.0f) const;
	virtual void glDraw(const RGBA& color) const {}

protected : 
	//////////////////////////////////////////////////////////////////////////
	// Basic Split Procedure
	class _SplitProc {};
	
	template<typename _QUADTREE>
	class SplitProc : public _SplitProc {
	public : 
		SplitProc(_QUADTREE *_qt) : qt(_qt) {}

		_CELL * createRootCell(unsigned int rootLevel) const;		

		inline const bool isSplit(const _CELL *cell) const;

		void splitCell(_CELL *cell) const;
		void postProcess() const;		

	public : 
		_QUADTREE *qt;
	};

public : 
	virtual _SplitProc* getSplitProc() const = 0;

public : 	
	unsigned int			maxLevel, rootLevel;
	_REAL					maxVal;

	Vector2<_REAL>			lc, uc, range;

	size_t					numCells;
	_CELL*					root;
	std::vector<_CELL*>		cells;
};

//////////////////////////////////////////////////////////////////////////
// Constructor/Destructor
template<typename _REAL, typename _CELL>
_QuadTree<_REAL, _CELL>::_QuadTree(unsigned int ml, Vector2<_REAL> lc, Vector2<_REAL> uc) : root(NULL) {	
	init(ml, lc, uc);	
}

template<typename _REAL, typename _CELL>
_QuadTree<_REAL, _CELL>::~_QuadTree() {
	release();
}

template<typename _REAL, typename _CELL>
void _QuadTree<_REAL, _CELL>::init() {
	init(maxLevel, lc, uc);
}

template<typename _REAL, typename _CELL>
void _QuadTree<_REAL, _CELL>::init(unsigned int _ml, Vector2<_REAL> _lc, Vector2<_REAL> _uc) {
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
void _QuadTree<_REAL, _CELL>::release() {
	cells.clear();

	if(root)
		geometryspace::release(root);	
	root = NULL;
}

//////////////////////////////////////////////////////////////////////////
// Utility
template<typename _REAL, typename _CELL>
void _QuadTree<_REAL, _CELL>::constrainedBdry(Vector2<_REAL> &pos) const {
	if ( pos[0] < lc[0] ) pos[0] = lc[0]; else if ( pos[0] > uc[0] ) pos[0] = uc[0];
	if ( pos[1] < lc[1] ) pos[1] = lc[1]; else if ( pos[1] > uc[1] ) pos[1] = uc[1];
}

//////////////////////////////////////////////////////////////////////////
// Get Location Code _ pos should lie between lc & uc
template<typename _REAL, typename _CELL>
inline Vector2<_REAL> _QuadTree<_REAL, _CELL>::getPos(const unsigned int loc[2]) const {	
	Vector2<_REAL> p(loc[0] / maxVal, loc[1] / maxVal);
	return (lc + range * p);							// does need Range Checker ?
}

template<typename _REAL, typename _CELL>
inline void _QuadTree<_REAL, _CELL>::getLocCode(unsigned int loc[2], const Vector2<_REAL> &pos) const {
	Vector2<_REAL> p = (pos - lc) / range;

	for(size_t i=0;i<2;++i) {
		loc[i] = (unsigned int)(p[i] * maxVal);
		if(loc[i] & (1 << rootLevel))
			--loc[i];
	}	
}

template<typename _REAL, typename _CELL>
inline Vector2<_REAL> _QuadTree<_REAL, _CELL>::getCellSize(unsigned int level) const {
	unsigned int divide = (1 << (rootLevel - level));
	return range / static_cast<_REAL>(divide);
}

template<typename _REAL, typename _CELL>
inline Vector2<_REAL> _QuadTree<_REAL, _CELL>::getCenterPos(const _CELL *cell) const {
	return getPos(cell->locCode) + getCellSize(cell->level) * static_cast<_REAL>(0.5f);
}

//////////////////////////////////////////////////////////////////////////
// Locate Function _ pos, p0, p1 should lie in the domain
template<typename _REAL, typename _CELL>
inline const _CELL * _QuadTree<_REAL, _CELL>::getLocateCell(const Vector2<_REAL> &pos) const {
	unsigned int loc[2];
	getLocCode(loc, pos);
	return QuadTreeLibrary::getLocateCell(root, loc);
}

template<typename _REAL, typename _CELL>
const _CELL * _QuadTree<_REAL, _CELL>::getLocateRegion(const Vector2<_REAL> &p0, const Vector2<_REAL> &p1) const {
	unsigned int loc0[2], loc1[2];
	getLocCode(loc0, p0);
	getLocCode(loc1, p1);
	
	unsigned int diff[2] = { loc0[0] ^ loc1[0], loc0[1] ^ loc1[1] };
	
	int level = rootLevel;
	int minLevel = rootLevel;
	while ((level != -1) && !(diff[0] & (1 << level))) level--;
	while ((minLevel > level) && !(diff[1] & (1 << minLevel))) minLevel--;
	
	return QuadTreeLibrary::getLocateCell(root, loc0, minLevel + 1);
}	

//////////////////////////////////////////////////////////////////////////
// Make Tree
template<typename _REAL, typename _CELL> template<typename SPLITPROC>
void _QuadTree<_REAL, _CELL>::makeTree(const SPLITPROC &sproc) {
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
void _QuadTree<_REAL, _CELL>::splitCell(_CELL *cell, const SPLITPROC &sproc) {
	if(!sproc.isSplit(cell))
		return;

	sproc.splitCell(cell);
}

//////////////////////////////////////////////////////////////////////////
// glDraw
template<typename _REAL, typename _CELL> 
void _QuadTree<_REAL, _CELL>::glDrawQuadTree(const RGBA &color, GLfloat width) const{
	for(unsigned int i=0;i<cells.size();++i) {
		Vector2<_REAL> pos = getPos(cells[i]->locCode);
		Vector2<_REAL> cellSize = getCellSize(cells[i]->level);

		glColor4f(color[0], color[1], color[2], color[3]);

		glLineWidth(width);
		glBegin(GL_LINE_STRIP);
			glVertex3f(pos[0]			 , pos[1]			 , 0.0f);
			glVertex3f(pos[0]+cellSize[0], pos[1]			 , 0.0f);
			glVertex3f(pos[0]+cellSize[0], pos[1]+cellSize[1], 0.0f);
			glVertex3f(pos[0]			 , pos[1]+cellSize[1], 0.0f);
			glVertex3f(pos[0]			 , pos[1]			 , 0.0f);
		glEnd();
	}	
}

//////////////////////////////////////////////////////////////////////////
// Basic Split Procedure
template<typename _REAL, typename _CELL> template<typename _QUADTREE>
_CELL * _QuadTree<_REAL, _CELL>::SplitProc<_QUADTREE>::createRootCell(unsigned int rootLevel) const {
	unsigned int rootCode[2] = {0,0};
	return new _CELL(rootCode, rootLevel, NULL);	
}

template<typename _REAL, typename _CELL> template<typename _QUADTREE>
inline const bool _QuadTree<_REAL, _CELL>::SplitProc<_QUADTREE>::isSplit(const _CELL *cell) const {
	return !(cell->level == 0);
}

template<typename _REAL, typename _CELL> template<typename _QUADTREE>
void _QuadTree<_REAL, _CELL>::SplitProc<_QUADTREE>::splitCell(_CELL *cell) const { 
	cell->children = new _CELL[4];
	for(unsigned int i=0;i<4;++i) {
		cell->children[i].parent = cell;
		cell->children[i].children = NULL;
		cell->children[i].level = cell->level - 1;

		for(int axis=0;axis<2;++axis) 
			cell->children[i].locCode[axis] = (cell->locCode[axis] | (((i & (1 << axis)) >> axis) << cell->children[i].level));

		qt->cells.push_back(&cell->children[i]);
	}
}

template<typename _REAL, typename _CELL> template<typename _QUADTREE>
void _QuadTree<_REAL, _CELL>::SplitProc<_QUADTREE>::postProcess() const {
	qt->numCells = qt->cells.size();
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Base QuadTree /////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
template<typename _REAL>
class QuadTree : public _QuadTree<_REAL, QTCell> {
public : 
	QuadTree(unsigned int ml = 2, Vector2<_REAL> lc = Vector2<_REAL>(0), Vector2<_REAL> uc = Vector2<_REAL>(1))
	: _QuadTree(ml, lc, uc) {}

public : 
	typename _QuadTree<_REAL, QTCell>::SplitProc<QuadTree> * getSplitProc() const;	
};

template<typename _REAL>
typename _QuadTree<_REAL, QTCell>::SplitProc<QuadTree<_REAL>> * QuadTree<_REAL>::getSplitProc() const {
	return new _QuadTree<_REAL, QTCell>::SplitProc<QuadTree<_REAL>>(const_cast<QuadTree<_REAL>*>(this));
}

};

#endif