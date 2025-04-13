#ifndef _GEOMETRY_QUADTREE_VERTICES_H_
#define _GEOMETRY_QUADTREE_VERTICES_H_

#include <hash_map>
#include "QuadTree.h"

namespace geometryspace {

//////////////////////////////////////////////////////////////////////////
// QuadTree_Vertices Base Class
template<typename _REAL, typename _CELL>
class _QuadTree_Vertices : public _QuadTree<_REAL, _CELL> {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	_QuadTree_Vertices(unsigned int ml, Vector2<_REAL> lc, Vector2<_REAL> uc);

	void init();
	void init(unsigned int ml, Vector2<_REAL> lc, Vector2<_REAL> uc);
	void release();

	const unsigned int vertIdx(unsigned int i, unsigned int j) const { return i + j * nVerticesAlongAxis;					}
	const unsigned int vertIdx(unsigned int locCode[2]) const		 { return locCode[0] + locCode[1] * nVerticesAlongAxis;	}

protected : 
	//////////////////////////////////////////////////////////////////////////
	// Basic Split Procedure
	template<typename _QUADTREE>
	class SplitProc : public _QuadTree<_REAL, _CELL>::SplitProc<_QUADTREE> {
	public : 
		SplitProc(_QUADTREE *qt) : _QuadTree<_REAL, _CELL>::SplitProc<_QUADTREE>(qt) {}

		_CELL * createRootCell(unsigned int rootLevel) const;

		inline const bool isSplit(_CELL *cell) const;

		void splitCell(_CELL *cell) const;
		void postProcess() const;
	};

public : 	
	size_t											numVertices;
	std::vector<Vector2<_REAL>>						posVerts;

protected :
	unsigned int									nVerticesAlongAxis;
	stdext::hash_map<unsigned int, unsigned int>	vertMap;

private : 
	//////////////////////////////////////////////////////////////////////////
	// variable for split routine
	unsigned int _vertUpdating[4][3];
};

//////////////////////////////////////////////////////////////////////////
// Constructor/Destructor
template<typename _REAL, typename _CELL>
_QuadTree_Vertices<_REAL, _CELL>::_QuadTree_Vertices(unsigned int _ml, Vector2<_REAL> _lc, Vector2<_REAL> _uc) 
: _QuadTree<_REAL, _CELL>(_ml, _lc, _uc) {
	init(_ml, _lc, _uc);

	// variable for split routine
	for(int i=0;i<4;++i) {
		int idx = 0;
		for(int j=0;j<4;++j) {
			if(i != (3 - j)) {
				_vertUpdating[i][idx++] = j;
			}
		}
	}
}

template<typename _REAL, typename _CELL>
void _QuadTree_Vertices<_REAL, _CELL>::init() {
	init(maxLevel, lc, uc);
}

template<typename _REAL, typename _CELL>
void _QuadTree_Vertices<_REAL, _CELL>::init(unsigned int _ml, Vector2<_REAL> _lc , Vector2<_REAL> _uc) {
	release();
	_QuadTree<_REAL, _CELL>::init(_ml, _lc, _uc);
	nVerticesAlongAxis = (1 << rootLevel) + 1;
}

template<typename _REAL, typename _CELL>
void _QuadTree_Vertices<_REAL, _CELL>::release() {
	_QuadTree<_REAL, _CELL>::release();
	posVerts.clear();
	vertMap.clear();	
}

//////////////////////////////////////////////////////////////////////////
// Basic Split Procedure
template<typename _REAL, typename _CELL> template<typename _QUADTREE>
_CELL * _QuadTree_Vertices<_REAL, _CELL>::SplitProc<_QUADTREE>::createRootCell(unsigned int rootLevel) const {
	unsigned int verts[4] = {0,1,2,3};	
	unsigned int locCodes[4][2] = {{0,0}, {1<<rootLevel,0}, {0,1<<rootLevel}, {1<<rootLevel,1<<rootLevel}};

	for(unsigned int i=0;i<4;++i) {
		qt->posVerts.push_back(qt->getPos(locCodes[i]));
		qt->vertMap[qt->vertIdx(locCodes[i])] = i;
	}

	return new _CELL(locCodes[0], rootLevel, NULL, verts);	
}

template<typename _REAL, typename _CELL> template<typename _QUADTREE>
inline const bool _QuadTree_Vertices<_REAL, _CELL>::SplitProc<_QUADTREE>::isSplit(_CELL *cell) const {
	return !(cell->level == 0);
}

template<typename _REAL, typename _CELL> template<typename _QUADTREE>
void _QuadTree_Vertices<_REAL, _CELL>::SplitProc<_QUADTREE>::splitCell(_CELL *cell) const { 
	_QuadTree<_REAL, _CELL>::SplitProc<_QUADTREE>::splitCell(cell);

	unsigned int binaryCellSize = (1 << (cell->children[0].level));

	Vector2<unsigned int> centerLocCode(cell->locCode[0] + binaryCellSize, cell->locCode[1] + binaryCellSize);
	Vector2<_REAL> centerPos(qt->getPos(centerLocCode.data));
	unsigned int vertIdx(qt->vertIdx(centerLocCode.data));

	unsigned int centerVert = static_cast<unsigned int>(qt->posVerts.size());
	qt->vertMap[vertIdx] = centerVert;
	qt->posVerts.push_back(centerPos);	

	unsigned int vertLocCode[2];
	for(size_t i=0;i<4;++i) {
		cell->children[i].vertices[3-i] = centerVert;

		for(size_t j=0;j<3;++j) {
			unsigned int vert = qt->_vertUpdating[i][j];

			for(size_t axis=0;axis<2;++axis) 
				vertLocCode[axis] = cell->children[i].locCode[axis] + (((vert & (1<<axis)) >> axis) << cell->children[0].level);

			vertIdx = qt->vertIdx(vertLocCode);
			if(!qt->vertMap.count(vertIdx)) {
				qt->vertMap[vertIdx] = cell->children[i].vertices[vert] = static_cast<unsigned int>(qt->posVerts.size());
				qt->posVerts.push_back(qt->getPos(vertLocCode));				
			}
			else cell->children[i].vertices[vert] = qt->vertMap[vertIdx];
		}
	}
}

template<typename _REAL, typename _CELL> template<typename _QUADTREE>
void _QuadTree_Vertices<_REAL, _CELL>::SplitProc<_QUADTREE>::postProcess() const {
	_QuadTree<_REAL, _CELL>::SplitProc<_QUADTREE>::postProcess();
	qt->numVertices = qt->posVerts.size();
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// QuadTree_Vertices /////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
template <typename _REAL>
class QuadTree_Vertices : public _QuadTree_Vertices<_REAL, QTCell_Vertices> {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	QuadTree_Vertices(unsigned int ml = 2, Vector2<_REAL> lc = Vector2<_REAL>(0), Vector2<_REAL> uc = Vector2<_REAL>(1))
	: _QuadTree_Vertices(ml, lc, uc) {}

public : 
	typename _QuadTree_Vertices<_REAL, QTCell_Vertices>::SplitProc<QuadTree_Vertices> * getSplitProc() const;
};

template<typename _REAL>
typename _QuadTree_Vertices<_REAL, QTCell_Vertices>::SplitProc<QuadTree_Vertices<_REAL>> * QuadTree_Vertices<_REAL>::getSplitProc() const {
	return new _QuadTree_Vertices<_REAL, QTCell_Vertices>::SplitProc<QuadTree_Vertices<_REAL>>(const_cast<QuadTree_Vertices<_REAL>*>(this));
}

};		// geometryspace

#endif
