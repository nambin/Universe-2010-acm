#ifndef _GEOMETRYSPACE_PLYLOADER_H_
#define _GEOMETRYSPACE_PLYLOADER_H_

#include <stdio.h>
#include <math.h>
#include <iostream>
#include "ply.h"
#include "../../../Basis/Basis/BasicLibrary.h"
#include "../../Principle.h"

namespace plyloader {

struct Vertex {
	float x,y,z;
	unsigned char r,g,b;
	float nx,ny,nz;	
	float confidence;
	//void *other_props;       // other properties 

	Vertex() : x(0), y(0), z(0), r(1), g(1), b(1), nx(0), ny(0), nz(0), confidence(0) {}
}; 

struct Face {
	unsigned char nverts;		// number of vertex indices in list 
	int *verts;					// vertex index list 
	//void *other_props;		// other properties 

	Face() : nverts(0), verts(NULL) {}
};

static char *elem_names[] = {	// list of the elements in the object 
	"vertex", "face"
};

static PlyProperty vert_props[] = {		// list of property information for a vertex 
	{"x", Float32, Float32, offsetof(Vertex,x), 0, 0, 0, 0},
	{"y", Float32, Float32, offsetof(Vertex,y), 0, 0, 0, 0},
	{"z", Float32, Float32, offsetof(Vertex,z), 0, 0, 0, 0},
	{"red", Uint8, Uint8, offsetof(Vertex,r), 0, 0, 0, 0},
	{"green", Uint8, Uint8, offsetof(Vertex,g), 0, 0, 0, 0},
	{"blue", Uint8, Uint8, offsetof(Vertex,b), 0, 0, 0, 0},	
	{"nx", Float32, Float32, offsetof(Vertex,nx), 0, 0, 0, 0},
	{"ny", Float32, Float32, offsetof(Vertex,ny), 0, 0, 0, 0},
	{"nz", Float32, Float32, offsetof(Vertex,nz), 0, 0, 0, 0},
	{"confidence", Float32, Float32, offsetof(Vertex,confidence), 0, 0, 0, 0},
};

static PlyProperty face_props[] = {		// list of property information for a face 
	{"vertex_indices", Int32, Int32, offsetof(Face,verts), 1, Uint8, Uint8, offsetof(Face,nverts)},
};

PRINCIPLEDLL std::ostream &operator<<(std::ostream &os, const Vertex &v);
PRINCIPLEDLL std::ostream &operator<<(std::ostream &os, const Face &v);
PRINCIPLEDLL void plyLoader(FILE *fp, int &nverts, int &nfaces, Vertex **&vlist, Face **&flist, int &per_vertex_color, int &has_normals);
PRINCIPLEDLL void plyConvert(const char *infile, const char *outfile);

};

#endif