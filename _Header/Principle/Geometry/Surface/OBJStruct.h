#ifndef _GEOMETRYSPACE_OBJSTRUCT_H_
#define _GEOMETRYSPACE_OBJSTRUCT_H_

#include "../../../Basis/Basis/BasicLibrary.h"
#include "../../../Basis/StringLibrary/StringLibrary.h"
#include "../../../Basis/DataPath/DataPath.h"
#include "../../../Basis/Color/RGBA.h"
#include "../../../Basis/Image/Image.h"
#include "../../Arithmetic/Vector/Vector2.h"
#include "../../Arithmetic/Vector/Vector3.h"
#include "../../Principle.h"

namespace geometryspace {

using std::vector;
using std::string;
using arithmetic::Vector2;
using arithmetic::Vector3;
using basis::RGBA;

struct PRINCIPLEDLL OBJMaterial {
public : 
	OBJMaterial() { init(); }
	void init();
	void set_DefaultMaterial();

public : 
	string		name;

	RGBA		diffuse;           
	RGBA		ambient;
	RGBA		specular;
	RGBA		emissive;
	float		shininess;        

	int			illum;
	string		map_Kd;		
	int			map_Kd_texture;
};

PRINCIPLEDLL std::ostream& operator<<(std::ostream &os, const OBJMaterial &mat);

struct OBJTriangle {
	int vIndex[3];
	int nIndex[3];
	int tIndex[3];
	int fIndex;
};

class OBJStruct;
struct PRINCIPLEDLL OBJGroup {
public : 
	OBJGroup(std::string _name = std::string("")) { init(_name); }
	void init(std::string _name);

	void setMaterial(int material);
	void pushTriangle(int idx);

	void glDraw(const OBJStruct &objstruct, int mode) const;

public : 
	std::string		name;          
	int				numTriangles;
	vector<int>		triangles;

	vector<int>		rangeMaterials;
	vector<int>		materials;

	OBJGroup*		next;          
};

class PRINCIPLEDLL OBJStruct {
public : 
	static const int OBJ_NONE		= 0;			// render with only vertices 
	static const int OBJ_FLAT		= (1 << 0);		// render with facet normals 
	static const int OBJ_SMOOTH		= (1 << 1);		// render with vertex normals 
	static const int OBJ_TEXTURE	= (1 << 2);     // render with texture coords 
	static const int OBJ_COLOR		= (1 << 3);		// render with colors 
	static const int OBJ_MATERIAL	= (1 << 4);		// render with materials 

public : 
	OBJStruct() { init(); }

	void init();
	void loadFromOBJ(const char* const filename);

	void setDrawMode(int mode) const;
	void setFaceNormals();

	void unitize();
	void translate(Vector3<float> trans);
	void scale(float s);

	void glDraw() const;

private : 
	void readMTL(const char* const name);
	OBJGroup* addGroup(const char* const name);
	OBJGroup* findGroup(const char* const name);
	const int findMaterial(const char* const name);
	const unsigned int makeTexture(const char* const name);

	void loadPass(const char* const name);

public : 
	std::string filename;
	std::string mtllibname;

	int numVertices, numTriangles;
	int numNormals, numTexCoords, numFaceNormals;
	int numMaterials, numGroups;

	vector<Vector3<float>>	posVerts;
	vector<OBJTriangle>		triangles;	

	vector<Vector3<float>>	normals;
	vector<Vector2<float>>	texCoords;
	vector<Vector3<float>>	faceNormals;

	vector<OBJMaterial>		materials;
	OBJGroup*				groups;

	Vector3<float>			position;

	mutable int				drawMode;
};

};		// geometryspace

#endif