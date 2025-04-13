#ifndef _GEOMETRYSPACE_PLYSTRUCT_H_
#define _GEOMETRYSPACE_PLYSTRUCT_H_

#include <vector>
#include "../../Arithmetic/Vector/Vector3.h"
#include "../../../Basis/Color/RGBA.h"
#include "../../Principle.h"

namespace geometryspace {

using namespace std;
using basis::RGBA;
using arithmetic::Vector3;

struct PLYVertex {
	Vector3<float> pos;
	Vector3<unsigned char> rgb;
	Vector3<float> normal;	
	float confidence;
};

struct PLYTriangle {
	int vIndex[3];
};

PRINCIPLEDLL std::ostream &operator<<(std::ostream &os, const PLYVertex &v);
PRINCIPLEDLL std::ostream &operator<<(std::ostream &os, const PLYTriangle &f);

class PRINCIPLEDLL PLYStruct {
public : 
	PLYStruct() { init(); }

	void init();
	void loadFromPLY(const char* const filename);

	void glDrawSimpleMesh(const RGBA &color = RGBA::white()) const;	
	void glDrawVertexColoredMesh() const;	

public : 
	int numVertices, numTriangles;

	vector<PLYVertex>		vertices;
	vector<PLYTriangle>		triangles;

	bool hasVertexColor, hasNormals;
};

};		// geometryspace

#endif