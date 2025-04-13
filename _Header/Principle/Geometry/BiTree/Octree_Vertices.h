#ifndef _GEOMETRY_OCTREE_VERTICES_H_
#define _GEOMETRY_OCTREE_VERTICES_H_

#include <hash_map>
#include "Octree.h"

namespace geometryspace {

//////////////////////////////////////////////////////////////////////////
// Octree_Vertices Base Class
template<typename _REAL, typename _CELL>
class _Octree_Vertices : public _Octree<_REAL, _CELL> {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	_Octree_Vertices(unsigned int ml, Vector3<_REAL> lc, Vector3<_REAL> uc);

	void init();
	void init(unsigned int ml, Vector3<_REAL> lc, Vector3<_REAL> uc);
	void release();

	const unsigned int vertIdx(unsigned int i, unsigned int j, unsigned int k) const { return i + (j + k * nVerticesAlongAxis) * nVerticesAlongAxis;							}
	const unsigned int vertIdx(unsigned int locCode[3]) const						 { return locCode[0] + (locCode[1] + locCode[2] * nVerticesAlongAxis) * nVerticesAlongAxis;	}

protected : 
	//////////////////////////////////////////////////////////////////////////
	// Basic Split Procedure
	template<typename _OCTREE>
	class SplitProc : public _Octree<_REAL, _CELL>::SplitProc<_OCTREE> {
	public : 
		SplitProc(_OCTREE *oct) : _Octree<_REAL, _CELL>::SplitProc<_OCTREE>(oct) {}

		_CELL * createRootCell(unsigned int rootLevel) const;

		inline const bool isSplit(_CELL *cell) const;

		void splitCell(_CELL *cell) const;
		void postProcess() const;
	};

public : 	
	size_t											numVertices;
	std::vector<Vector3<_REAL>>						posVerts;

protected :
	unsigned int									nVerticesAlongAxis;
	stdext::hash_map<unsigned int, unsigned int>	vertMap;

private : 
	//////////////////////////////////////////////////////////////////////////
	// variable for split routine
	unsigned int _vertUpdating[8][7];
};

//////////////////////////////////////////////////////////////////////////
// Constructor/Destructor
template<typename _REAL, typename _CELL>
_Octree_Vertices<_REAL, _CELL>::_Octree_Vertices(unsigned int ml, Vector3<_REAL> lc, Vector3<_REAL> uc) 
: _Octree<_REAL, _CELL>(ml, lc, uc) 
{
	init(ml, lc, uc);

	// variable for split routine
	for(int i=0;i<8;++i) {
		int idx = 0;
		for(int j=0;j<8;++j) {
			if(i != (7 - j)) {
				_vertUpdating[i][idx++] = j;
			}
		}
	}
}

template<typename _REAL, typename _CELL>
void _Octree_Vertices<_REAL, _CELL>::init() {
	init(maxLevel, lc, uc);
}

template<typename _REAL, typename _CELL>
void _Octree_Vertices<_REAL, _CELL>::init(unsigned int ml, Vector3<_REAL> lc , Vector3<_REAL> uc) {
	release();
	_Octree<_REAL, _CELL>::init(ml, lc, uc);
	nVerticesAlongAxis = (1 << rootLevel) + 1;
}

template<typename _REAL, typename _CELL>
void _Octree_Vertices<_REAL, _CELL>::release() {
	_Octree<_REAL, _CELL>::release();
	posVerts.clear();
	vertMap.clear();	
}

//////////////////////////////////////////////////////////////////////////
// Basic Split Procedure
template<typename _REAL, typename _CELL> template<typename _QUADTREE>
_CELL * _Octree_Vertices<_REAL, _CELL>::SplitProc<_QUADTREE>::createRootCell(unsigned int rootLevel) const {
	unsigned int verts[8] = {0,1,2,3,4,5,6,7};	
	unsigned int locCodes[8][3] = {{0,0,0           }, {1<<rootLevel,0,0           }, {0,1<<rootLevel,0           }, {1<<rootLevel,1<<rootLevel,0           }
								 , {0,0,1<<rootLevel}, {1<<rootLevel,0,1<<rootLevel}, {0,1<<rootLevel,1<<rootLevel}, {1<<rootLevel,1<<rootLevel,1<<rootLevel}};

	for(unsigned int i=0;i<8;++i) {
		oct->posVerts.push_back(oct->getPos(locCodes[i]));
		oct->vertMap[oct->vertIdx(locCodes[i])] = i;
	}

	return new _CELL(locCodes[0], rootLevel, NULL, verts);	
}

template<typename _REAL, typename _CELL> template<typename _QUADTREE>
inline const bool _Octree_Vertices<_REAL, _CELL>::SplitProc<_QUADTREE>::isSplit(_CELL *cell) const {
	//return !(cell->level == 0);
	if(cell->level == 0)
		return false;
	return (rand() / (double)RAND_MAX) < 0.4f;
}

template<typename _REAL, typename _CELL> template<typename _QUADTREE>
void _Octree_Vertices<_REAL, _CELL>::SplitProc<_QUADTREE>::splitCell(_CELL *cell) const { 
	_Octree<_REAL, _CELL>::SplitProc<_QUADTREE>::splitCell(cell);

	unsigned int binaryCellSize = (1 << (cell->children[0].level));

	Vector3<unsigned int> centerLocCode(cell->locCode[0] + binaryCellSize, cell->locCode[1] + binaryCellSize, cell->locCode[2] + binaryCellSize);
	Vector3<_REAL> centerPos(oct->getPos(centerLocCode.data));
	unsigned int vertIdx(oct->vertIdx(centerLocCode.data));

	unsigned int centerVert = static_cast<unsigned int>(oct->posVerts.size());
	oct->vertMap[vertIdx] = centerVert;
	oct->posVerts.push_back(centerPos);	

	unsigned int vertLocCode[3];
	for(size_t i=0;i<8;++i) {
		cell->children[i].vertices[7-i] = centerVert;

		for(size_t j=0;j<7;++j) {
			unsigned int vert = oct->_vertUpdating[i][j];

			for(unsigned int axis=0;axis<3;++axis) 
				vertLocCode[axis] = cell->children[i].locCode[axis] + (((vert & (1<<axis)) >> axis) << cell->children[0].level);

			vertIdx = oct->vertIdx(vertLocCode);
			if(!oct->vertMap.count(vertIdx)) {
				oct->vertMap[vertIdx] = cell->children[i].vertices[vert] = static_cast<unsigned int>(oct->posVerts.size());
				oct->posVerts.push_back(oct->getPos(vertLocCode));				
			}
			else cell->children[i].vertices[vert] = oct->vertMap[vertIdx];
		}
	}
}

template<typename _REAL, typename _CELL> template<typename _QUADTREE>
void _Octree_Vertices<_REAL, _CELL>::SplitProc<_QUADTREE>::postProcess() const {
	_Octree<_REAL, _CELL>::SplitProc<_QUADTREE>::postProcess();
	oct->numVertices = oct->posVerts.size();
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Octree_Vertices /////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
template <typename _REAL>
class Octree_Vertices : public _Octree_Vertices<_REAL, OctCell_Vertices> {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Octree_Vertices(unsigned int ml = 2, Vector3<_REAL> lc = Vector3<_REAL>(0), Vector3<_REAL> uc = Vector3<_REAL>(1))
	: _Octree_Vertices(ml, lc, uc) {}

public : 
	typename _Octree_Vertices<_REAL, OctCell_Vertices>::SplitProc<Octree_Vertices> * getSplitProc() const;
};

template<typename _REAL>
typename _Octree_Vertices<_REAL, OctCell_Vertices>::SplitProc<Octree_Vertices<_REAL>> * Octree_Vertices<_REAL>::getSplitProc() const {
	return new _Octree_Vertices<_REAL, OctCell_Vertices>::SplitProc<Octree_Vertices<_REAL>>(const_cast<Octree_Vertices<_REAL>*>(this));
}

};		// geometryspace

#endif
