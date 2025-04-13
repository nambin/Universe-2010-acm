#ifndef _GEOMETRY_QUADTREE_CELL_LIBRARY_H_
#define _GEOMETRY_QUADTREE_CELL_LIBRARY_H_

#include <cmath>
#include <vector>
#include "../../Arithmetic/Vector/Vector2.h"
#include "../../VirtualSpace/RenderingPipeline/RenderingPipelineSetting.h"

namespace geometryspace {

using basis::RGBA;
using arithmetic::Vector2;

class QuadTreeLibrary {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Cell property
	template<typename _QTCELL>
	inline static const bool isLeaf(const _QTCELL *cell) { 
		return !cell->children; 
	}

	//////////////////////////////////////////////////////////////////////////
	// Release QuadTree Cell
	template<typename _QTCELL>
	static void release(_QTCELL* cell) {
		if(!cell->children)
			return;

		for(unsigned int i=0;i<4;++i) 
			release(&cell->children[i]);
		delete[] cell->children;
	}

	//////////////////////////////////////////////////////////////////////////
	// getLocateCell - loc[2] should lie in this cell
	template<typename _QTCELL>
	static const _QTCELL * getLocateCell(const _QTCELL *_cell, const unsigned int loc[2], unsigned int destLevel = 0) {
		if(isLeaf(_cell) || destLevel >= _cell->level)
			return _cell;

		unsigned int nextLevel = _cell->level - 1;
		unsigned int nIter = nextLevel - destLevel + 1;
		const _QTCELL * cell = _cell;

		for(unsigned int i=0;i<nIter;++i) {
			unsigned int childBranchBit = 1 << (nextLevel); 
			unsigned int childIndex = ((loc[0] & childBranchBit) >> nextLevel) 
				+ (((loc[1] & childBranchBit) << 1) >> nextLevel--); 

			cell = &((cell->children)[childIndex]); 
			if(isLeaf(cell)) break; 
		} 

		return cell;
	}

	template<typename _QTCELL>
	static const _QTCELL * getCommonAncestor(const _QTCELL *_cell, const unsigned int loc, unsigned int axis) {
		unsigned int binaryDiff = _cell->locCode[axis] ^ loc;
		unsigned int cellLevel = _cell->level;
		const _QTCELL * cell = _cell;

		while (binaryDiff & (1 << cellLevel)) { 
			cell = cell->parent; 
			++cellLevel; 
		} 

		return cell;
	}

	//////////////////////////////////////////////////////////////////////////
	// Neighbor information
	template<typename _QTCELL>
	inline static const bool isXNeighbor(const _QTCELL *qc0, const _QTCELL *qc1) { 
		return (isRightNeighbor(qc0, qc1) || isRightNeighbor(qc1, qc0));	
	}

	template<typename _QTCELL>
	inline static const bool isYNeighbor(const _QTCELL *qc0, const _QTCELL *qc1) { 
		return (isUpNeighbor(qc0, qc1) || isUpNeighbor(qc1, qc0));			
	}

	template<typename _QTCELL>
	inline static const bool isRightNeighbor(const _QTCELL *cell, const _QTCELL *right) {
		if(right->locCode[0] == cell->locCode[0] + (1 << cell->level))
			return isRangeOverlap(cell, right, 1);
		return false;
	}

	template<typename _QTCELL>
	inline static const bool isUpNeighbor(const _QTCELL *cell, const _QTCELL *up) {
		if(up->locCode[1] == cell->locCode[1] + (1 << cell->level))
			return isRangeOverlap(cell, up, 0);
		return false;
	}

	template<typename _QTCELL>
	inline static const bool isRangeOverlap(const _QTCELL *qc0, const _QTCELL *qc1, unsigned int axis) {
		if(qc0->level > qc1->level) return (qc0->locCode[axis] <= qc1->locCode[axis] && qc1->locCode[axis] < (qc0->locCode[axis] + (1 << qc0->level)));
		else						return (qc1->locCode[axis] <= qc0->locCode[axis] && qc0->locCode[axis] < (qc1->locCode[axis] + (1 << qc1->level)));
	}

	template<typename _QTCELL>
	static void delete_NonLeafNeighbors(_QTCELL *cell) {
		if(!isLeaf(cell))
			cell->neighbors.clear();

		std::vector<_QTCELL*> _neighbor;
		for(unsigned int i=0;i<cell->neighbors.size();++i) {
			if(isLeaf(cell->neighbors[i])) {
				_neighbor.push_back(cell->neighbors[i]);
			}
		}
		cell->neighbors = _neighbor;
	}

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	template<typename _QTCELL, typename _QUADTREE>
	static void glDrawCell(const _QTCELL *cell, const _QUADTREE *qt, const RGBA &color, GLfloat width = 1.0f) {
		Vector2<float> pos = Vector2<float>(qt->getPos(cell->locCode));
		Vector2<float> cellSize = Vector2<float>(qt->getCellSize(cell->level));

		glColor4f(color[0], color[1], color[2], color[3]);

		glLineWidth(width);
		glBegin(GL_LINE_STRIP);
		glVertex3f(pos[0]			 , pos[1]			 , 0.001f);
		glVertex3f(pos[0]+cellSize[0], pos[1]			 , 0.001f);
		glVertex3f(pos[0]+cellSize[0], pos[1]+cellSize[1], 0.001f);
		glVertex3f(pos[0]			 , pos[1]+cellSize[1], 0.001f);
		glVertex3f(pos[0]			 , pos[1]			 , 0.001f);
		glEnd();
	}

	template<typename _QTCELL, typename _QUADTREE>
	static void glDrawCellVertices(const _QTCELL *cell, const _QUADTREE *qt, const RGBA &color, GLfloat size = 5.0f) {

		glColor4f(color[0], color[1], color[2], color[3]);

		glPointSize(size);
		glBegin(GL_POINTS);
		for(unsigned int i=0;i<4;++i) {
			Vector2<float> pos = Vector2<float>(qt->posVerts[cell->vertices[i]]);
			glVertex3f(pos[0], pos[1], 0.001f);
		}
		glEnd();
	}

	template<typename _QTCELL, typename _QUADTREE>
	static void glDrawCellNeighbors(const _QTCELL *cell, const _QUADTREE *qt, const RGBA &color, GLfloat width = 1.0f) {

		Vector2<float> cellCenter = qt->getCenterPos(cell);

		glColor4f(color[0], color[1], color[2], color[3]);

		glLineWidth(width);
		glBegin(GL_LINES);
		for(unsigned int i=0;i<cell->neighbors.size();++i) {
			Vector2<float> neighborCenter = qt->getCenterPos(cell->neighbors[i]);

			glVertex3f(cellCenter[0], cellCenter[1], 0.001f);
			glVertex3f(neighborCenter[0], neighborCenter[1], 0.001f);
		}
		glEnd();
	}
};

};		// geometryspace

#endif