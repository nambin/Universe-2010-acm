#ifndef _GEOMETRY_QUADTREE_NEIGHBORS_H_
#define _GEOMETRY_QUADTREE_NEIGHBORS_H_

#include <hash_map>
#include "QuadTree_Vertices.h"

namespace geometryspace {

//////////////////////////////////////////////////////////////////////////
// QuadTree_Neighbors Base Class
template<typename _REAL, typename _CELL>
class _QuadTree_Neighbors : public _QuadTree_Vertices<_REAL, _CELL> {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	_QuadTree_Neighbors(unsigned int ml, Vector2<_REAL> lc, Vector2<_REAL> uc);

protected : 
	//////////////////////////////////////////////////////////////////////////
	// Basic Split Procedure
	template<typename _QUADTREE>
	class SplitProc : public _QuadTree_Vertices<_REAL, _CELL>::SplitProc<_QUADTREE> {
	public : 
		SplitProc(_QUADTREE *qt) : _QuadTree_Vertices<_REAL, _CELL>::SplitProc<_QUADTREE>(qt) {}

		_CELL * createRootCell(unsigned int rootLevel) const;
		
		inline const bool isSplit(_CELL *cell) const;
		
		void splitCell(_CELL *cell) const;
		void postProcess() const;
	};
};

//////////////////////////////////////////////////////////////////////////
// Constructor/Destructor
template<typename _REAL, typename _CELL>
_QuadTree_Neighbors<_REAL,_CELL>::_QuadTree_Neighbors(unsigned int ml, Vector2<_REAL> lc, Vector2<_REAL> uc) 
: _QuadTree_Vertices<_REAL,_CELL>(ml, lc, uc) {}

//////////////////////////////////////////////////////////////////////////
// Basic Split Procedure
template<typename _REAL, typename _CELL> template<typename _QUADTREE>
_CELL * _QuadTree_Neighbors<_REAL, _CELL>::SplitProc<_QUADTREE>::createRootCell(unsigned int rootLevel) const {
	return _QuadTree_Vertices<_REAL, _CELL>::SplitProc<_QUADTREE>::createRootCell(rootLevel);	
}

template<typename _REAL, typename _CELL> template<typename _QUADTREE>
inline const bool _QuadTree_Neighbors<_REAL, _CELL>::SplitProc<_QUADTREE>::isSplit(_CELL *cell) const {
	return !(cell->level == 0);
}

template<typename _REAL, typename _CELL> template<typename _QUADTREE>
void _QuadTree_Neighbors<_REAL, _CELL>::SplitProc<_QUADTREE>::splitCell(_CELL *cell) const { 
	_QuadTree_Vertices<_REAL, _CELL>::SplitProc<_QUADTREE>::splitCell(cell);

	std::vector<_CELL*> &neighbors = cell->neighbors;

	for(size_t i=0;i<4;++i) {
		for(size_t axis=0;axis<2;++axis) {
			unsigned int idxNeighbor = (i ^ (1 << axis)) & 3;
			cell->children[i].neighbors.push_back(&cell->children[idxNeighbor]);
		}

		for(size_t idx=0;idx<neighbors.size();++idx) {
			if(QuadTreeLibrary::isXNeighbor(&cell->children[i], neighbors[idx]) 
			|| QuadTreeLibrary::isYNeighbor(&cell->children[i], neighbors[idx])) {
				cell->children[i].addNeighbor(neighbors[idx]);
				neighbors[idx]->addNeighbor(&cell->children[i]);
			}
		}
	}
}

template<typename _REAL, typename _CELL> template<typename _QUADTREE>
void _QuadTree_Neighbors<_REAL, _CELL>::SplitProc<_QUADTREE>::postProcess() const { 
	_QuadTree_Vertices<_REAL, _CELL>::SplitProc<_QUADTREE>::postProcess();
	for(size_t i=0;i<qt->cells.size();++i) {
		QuadTreeLibrary::delete_NonLeafNeighbors(qt->cells[i]);
	}
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// QuadTree_Neighbors ////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
template <typename _REAL>
class QuadTree_Neighbors : public _QuadTree_Neighbors<_REAL, QTCell_Neighbors> {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	QuadTree_Neighbors(unsigned int ml = 2, Vector2<_REAL> lc = Vector2<_REAL>(0), Vector2<_REAL> uc = Vector2<_REAL>(1))
	: _QuadTree_Neighbors(ml, lc, uc) {}

public : 
	typename _QuadTree_Neighbors<_REAL, QTCell_Neighbors>::SplitProc<QuadTree_Neighbors> * getSplitProc() const;
};

template<typename _REAL>
typename _QuadTree_Neighbors<_REAL, QTCell_Neighbors>::SplitProc<QuadTree_Neighbors<_REAL>> * QuadTree_Neighbors<_REAL>::getSplitProc() const {
	return new _QuadTree_Neighbors<_REAL, QTCell_Neighbors>::SplitProc<QuadTree_Neighbors<_REAL>>(const_cast<QuadTree_Neighbors<_REAL>*>(this));
}

};		// geometryspace

#endif