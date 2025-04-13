#ifndef _GEOMETRY_QUADTREE_CELL_H_
#define _GEOMETRY_QUADTREE_CELL_H_

#include "../../Principle.h"
#include "QuadTreeCell_Library.h"

namespace geometryspace {

class PRINCIPLEDLL QTCell {
public : 
	QTCell() : parent(NULL), children(NULL) {}	
	QTCell(unsigned int _locCode[2], unsigned int _level, QTCell* _parent) 
	: level(_level), parent(_parent), children(NULL) { 
		memcpy(locCode, _locCode, 2 * sizeof(unsigned int));
	}
	
public : 
	unsigned int	locCode[2];
	unsigned int	level;		
	QTCell*			parent;		
	QTCell*			children;		
};

class PRINCIPLEDLL QTCell_Vertices {
public : 
	QTCell_Vertices() : parent(NULL), children(NULL) {}	
	QTCell_Vertices(unsigned int _locCode[2], unsigned int _level, QTCell_Vertices* _parent, unsigned int vert[4])
	: level(_level), parent(_parent), children(NULL) { 
		memcpy(locCode, _locCode, 2 * sizeof(unsigned int));
		memcpy(vertices, vert, 4 * sizeof(unsigned int));
	}

public : 
	unsigned int		locCode[2];
	unsigned int		level;	
	QTCell_Vertices*	parent;	
	QTCell_Vertices*	children;

	unsigned int		vertices[4];
};

class PRINCIPLEDLL QTCell_Neighbors {
public : 
	QTCell_Neighbors() : parent(NULL), children(NULL) {}	
	QTCell_Neighbors(unsigned int _locCode[2], unsigned int _level, QTCell_Neighbors* _parent, unsigned int vert[4])
	: level(_level), parent(_parent), children(NULL) { 
		memcpy(locCode, _locCode, 2 * sizeof(unsigned int));
		memcpy(vertices, vert, 4 * sizeof(unsigned int));
	}

	void addNeighbor(QTCell_Neighbors *cell) { neighbors.push_back(cell); }

public : 
	unsigned int		locCode[2];
	unsigned int		level;	
	QTCell_Neighbors*	parent;	
	QTCell_Neighbors*	children;

	unsigned int		vertices[4];

	std::vector<QTCell_Neighbors*> neighbors;
};

};		// geometryspace

#endif