#include <stdafx.h>
#include <Principle/Geometry/Surface/OBJStruct.h>

namespace geometryspace {

void OBJMaterial::init() {
	name.clear();

	ambient = diffuse = specular = emissive = RGBA::black();	
	shininess = 0.0f;

	illum = 2;
	map_Kd.clear();
	map_Kd_texture = -1;
}

void OBJMaterial::set_DefaultMaterial() {
	name = string("default");

	diffuse = RGBA::white() * 0.8f;
	ambient = RGBA::white() * 0.4f;
	specular = emissive = RGBA::black();		
	shininess = 0.0f;

	illum = 2;
	map_Kd.clear();
	map_Kd_texture = -1;
}

std::ostream& operator<<(std::ostream &os, const OBJMaterial &mat) {
	os << "Name _ " << mat.name << std::endl;
	os << "Diffuse _ " << mat.diffuse << std::endl;
	os << "Ambient _ " << mat.ambient << std::endl;
	os << "Specular _ " << mat.specular << std::endl;
	os << "Emissive _ " << mat.emissive << std::endl;
	os << "Shininess _ " << mat.shininess << std::endl;
	os << "Type of Illumination _ " << mat.illum << std::endl;
	os << "Map Kd _ " << mat.map_Kd.c_str() << std::endl;
	return os;
}

void OBJGroup::init(std::string _name) {

	name = _name;
	numTriangles = 0;
	triangles.clear();

	rangeMaterials.clear(); rangeMaterials.push_back(0); rangeMaterials.push_back(0); 
	materials.clear();		materials.push_back(0);

	next = NULL;
}

void OBJGroup::setMaterial(int material) {

	if(materials[materials.size()-1] == material)
		return;

	if(rangeMaterials[materials.size()-1] == rangeMaterials[materials.size()]) {
		materials[materials.size()-1] = material;
		return;
	}
	
	rangeMaterials.push_back(rangeMaterials[materials.size()]);
	materials.push_back(material);
}

void OBJGroup::pushTriangle(int idx) {

	triangles.push_back(idx);
	++numTriangles;	
	++rangeMaterials[rangeMaterials.size()-1];
}

void OBJStruct::init() {

	filename.clear();
	mtllibname.clear();
	numVertices = numTriangles = numNormals = 0;
	numTexCoords = numFaceNormals = numMaterials = numGroups = 0;

	posVerts.clear();
	triangles.clear();
	normals.clear();
	texCoords.clear();
	faceNormals.clear();
	materials.clear();
	groups = NULL;
	position = Vector3<float>(static_cast<float>(0));

	drawMode = (OBJ_FLAT | OBJ_SMOOTH | OBJ_TEXTURE | OBJ_COLOR | OBJ_MATERIAL);
}

void OBJStruct::setDrawMode(int mode) const {
	drawMode = mode;
}

};		// geometryspace
