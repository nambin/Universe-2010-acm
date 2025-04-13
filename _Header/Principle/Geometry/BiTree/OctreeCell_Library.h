#ifndef _GEOMETRY_OCTREE_CELL_LIBRARY_H_
#define _GEOMETRY_OCTREE_CELL_LIBRARY_H_

#include <cmath>
#include <vector>
#include "../../Arithmetic/Vector/Vector3.h"
#include "../../VirtualSpace/RenderingPipeline/RenderingPipelineSetting.h"

namespace geometryspace {

using basis::RGBA;
using arithmetic::Vector3;

class OctreeLibrary {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Cell property
	template<typename _OCTCELL>
	inline static const bool isLeaf(const _OCTCELL *cell) { 
		return !cell->children; 
	}

	//////////////////////////////////////////////////////////////////////////
	// Release QuadTree Cell
	template<typename _OCTCELL>
	static void release(_OCTCELL* cell) {
		if(!cell->children)
			return;

		for(size_t i=0;i<8;++i) 
			release(&cell->children[i]);
		delete[] cell->children;
	}

	//////////////////////////////////////////////////////////////////////////
	// getLocateCell - loc[3] should lie in this cell
	template<typename _OCTCELL>
	static const _OCTCELL * getLocateCell(const _OCTCELL *_cell, const unsigned int loc[3], unsigned int destLevel = 0) {
		if(isLeaf(_cell) || destLevel >= _cell->level)
			return _cell;

		unsigned int nextLevel = _cell->level - 1;
		unsigned int nIter = nextLevel - destLevel + 1;
		const _OCTCELL * cell = _cell;

		for(unsigned int i=0;i<nIter;++i) {
			unsigned int childBranchBit = 1 << (nextLevel); 
			unsigned int childIndex = ((loc[0] & childBranchBit) >> nextLevel) 
				+ (((loc[1] & childBranchBit) << 1) >> nextLevel)
				+ (((loc[2] & childBranchBit) << 2) >> nextLevel--); 

			cell = &((cell->children)[childIndex]); 
			if(isLeaf(cell)) break; 
		} 

		return cell;
	}

	template<typename _OCTCELL>
	static const _OCTCELL * getCommonAncestor(const _OCTCELL *_cell, const unsigned int loc, unsigned int axis) {
		unsigned int binaryDiff = _cell->locCode[axis] ^ loc;
		unsigned int cellLevel = _cell->level;
		const _OCTCELL * cell = _cell;

		while (binaryDiff & (1 << cellLevel)) { 
			cell = cell->parent; 
			++cellLevel; 
		} 

		return cell;
	}

	//////////////////////////////////////////////////////////////////////////
	// Neighbor information
	template<typename _OCTCELL>
	inline static const bool isXNeighbor(const _OCTCELL *qc0, const _OCTCELL *qc1) { 
		return (isXPlusNeighbor(qc0, qc1) || isXPlusNeighbor(qc1, qc0));	
	}

	template<typename _OCTCELL>
	inline static const bool isYNeighbor(const _OCTCELL *qc0, const _OCTCELL *qc1) { 
		return (isYPlusNeighbor(qc0, qc1) || isYPlusNeighbor(qc1, qc0));			
	}

	template<typename _OCTCELL>
	inline static const bool isZNeighbor(const _OCTCELL *qc0, const _OCTCELL *qc1) { 
		return (isZPlusNeighbor(qc0, qc1) || isZPlusNeighbor(qc1, qc0));			
	}

	template<typename _OCTCELL>
	inline static const bool isXPlusNeighbor(const _OCTCELL *cell, const _OCTCELL *up) {
		if(right->locCode[0] == cell->locCode[0] + (1 << cell->level))
			return (isRangeOverlap(cell, up, 1) && isRangeOverlap(cell, up, 2));
		return false;
	}

	template<typename _OCTCELL>
	inline static const bool isYPlusNeighbor(const _OCTCELL *cell, const _OCTCELL *up) {
		if(up->locCode[1] == cell->locCode[1] + (1 << cell->level))
			return (isRangeOverlap(cell, up, 0) && isRangeOverlap(cell, up, 2));
		return false;
	}

	template<typename _OCTCELL>
	inline static const bool isZPlusNeighbor(const _OCTCELL *cell, const _OCTCELL *up) {
		if(up->locCode[2] == cell->locCode[2] + (1 << cell->level))
			return (isRangeOverlap(cell, up, 0) && isRangeOverlap(cell, up, 1));
		return false;
	}

	template<typename _OCTCELL>
	inline static const bool isRangeOverlap(const _OCTCELL *qc0, const _OCTCELL *qc1, unsigned int axis) {
		if(qc0->level > qc1->level) return (qc0->locCode[axis] <= qc1->locCode[axis] && qc1->locCode[axis] < (qc0->locCode[axis] + (1 << qc0->level)));
		else						return (qc1->locCode[axis] <= qc0->locCode[axis] && qc0->locCode[axis] < (qc1->locCode[axis] + (1 << qc1->level)));
	}

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	template<typename _OCTCELL, typename _OCTREE>
	static void glDrawCell(const _OCTCELL *cell, const _OCTREE *qt, const RGBA &color, GLfloat width = 1.0f) {
		Vector3<float> pos = Vector3<float>(qt->getPos(cell->locCode));
		Vector3<float> cellSize = Vector3<float>(qt->getCellSize(cell->level));

		glColor4f(color[0], color[1], color[2], color[3]);

		glLineWidth(width);
		glBegin(GL_LINE_STRIP);
		glVertex3f(pos[0]			 , pos[1]			 , pos[2]);
		glVertex3f(pos[0]+cellSize[0], pos[1]			 , pos[2]);
		glVertex3f(pos[0]+cellSize[0], pos[1]+cellSize[1], pos[2]);
		glVertex3f(pos[0]			 , pos[1]+cellSize[1], pos[2]);
		glVertex3f(pos[0]			 , pos[1]			 , pos[2]);
		glEnd();

		glBegin(GL_LINE_STRIP);
		glVertex3f(pos[0]			 , pos[1]			 , pos[2]+cellSize[2]);
		glVertex3f(pos[0]+cellSize[0], pos[1]			 , pos[2]+cellSize[2]);
		glVertex3f(pos[0]+cellSize[0], pos[1]+cellSize[1], pos[2]+cellSize[2]);
		glVertex3f(pos[0]			 , pos[1]+cellSize[1], pos[2]+cellSize[2]);
		glVertex3f(pos[0]			 , pos[1]			 , pos[2]+cellSize[2]);
		glEnd();

		glBegin(GL_LINES);
		glVertex3f(pos[0]			 , pos[1]			 , pos[2]);
		glVertex3f(pos[0]			 , pos[1]			 , pos[2]+cellSize[2]);

		glVertex3f(pos[0]+cellSize[0], pos[1]			 , pos[2]);
		glVertex3f(pos[0]+cellSize[0], pos[1]			 , pos[2]+cellSize[2]);

		glVertex3f(pos[0]+cellSize[0], pos[1]+cellSize[1], pos[2]);
		glVertex3f(pos[0]+cellSize[0], pos[1]+cellSize[1], pos[2]+cellSize[2]);

		glVertex3f(pos[0]			 , pos[1]+cellSize[1], pos[2]);
		glVertex3f(pos[0]			 , pos[1]+cellSize[1], pos[2]+cellSize[2]);
		glEnd();
	}

	template<typename _OCTCELL, typename _OCTREE>
	static void glDrawCellVertices(const _OCTCELL *cell, const _OCTREE *qt, const RGBA &color, GLfloat size = 5.0f) {
		glColor4f(color[0], color[1], color[2], color[3]);

		glPointSize(size);
		glBegin(GL_POINTS);
		for(unsigned int i=0;i<8;++i) {
			Vector3<float> pos = Vector3<float>(qt->posVerts[cell->vertices[i]]);
			glVertex3f(pos[0], pos[1], pos[2]);
		}
		glEnd();
	}
};

};		// geometryspace

#endif