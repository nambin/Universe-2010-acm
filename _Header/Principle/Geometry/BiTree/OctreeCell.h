#ifndef _GEOMETRY_OCTREE_CELL_H_
#define _GEOMETRY_OCTREE_CELL_H_

#include "../../Principle.h"
#include "OctreeCell_Library.h"

namespace geometryspace {

class PRINCIPLEDLL OctCell {
public : 
	OctCell() : parent(NULL), children(NULL) {}	
	OctCell(unsigned int _locCode[3], unsigned int _level, OctCell* _parent) 
	: level(_level), parent(_parent), children(NULL) { 
		memcpy(locCode, _locCode, 3 * sizeof(unsigned int));
	}

public : 
	unsigned int	locCode[3];
	unsigned int	level;		
	OctCell*		parent;		
	OctCell*		children;		
};

class PRINCIPLEDLL OctCell_Vertices {
public : 
	OctCell_Vertices() : parent(NULL), children(NULL) {}	
	OctCell_Vertices(unsigned int _locCode[3], unsigned int _level, OctCell_Vertices* _parent, unsigned int vert[8])
	: level(_level), parent(_parent), children(NULL) { 
		memcpy(locCode, _locCode, 3 * sizeof(unsigned int));
		memcpy(vertices, vert, 8 * sizeof(unsigned int));
	}

public : 
	unsigned int		locCode[3];
	unsigned int		level;	
	OctCell_Vertices*	parent;	
	OctCell_Vertices*	children;

	unsigned int		vertices[8];
};

};		// geometryspace

#endif