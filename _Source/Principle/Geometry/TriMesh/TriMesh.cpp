#include <stdafx.h>
#include <Principle/Geometry/Surface/TriMesh.h>

using namespace geometryspace;

//////////////////////////////////////////////////////////////////////////
// MeshMaterial
void MeshMaterial::init() {
	material.init();
	map_Kd_filename.clear();
	map_Kd_texture = -1;
}

void MeshMaterial::set(const MeshMaterial &mat) {
	(*this) = mat;
}

void MeshMaterial::set(const OBJMaterial &mat) {
	material.setAmbient(mat.ambient);
	material.setDiffuse(mat.diffuse);
	material.setSpecular(mat.specular);
	material.setEmissive(mat.emissive);
	material.setShininess(mat.shininess);
	map_Kd_filename = mat.map_Kd;
	map_Kd_texture = mat.map_Kd_texture;
}

void MeshMaterial::set_DefaultMaterial() {
	OBJMaterial mat;
	mat.set_DefaultMaterial();
	set(mat);	
}

MeshMaterial& MeshMaterial::operator=(const OBJMaterial &mat) {
	set(mat);
	return (*this);
}

//////////////////////////////////////////////////////////////////////////
// MeshGroup
void MeshGroup::init(std::string _name) {
	name = _name;
	numTriangles = 0;
	triangles.clear();

	rangeMaterials.clear();
	materials.clear();
}

void MeshGroup::set(const MeshGroup &group) {
	(*this) = group;
}

void MeshGroup::set(const OBJGroup &group) {
	name = group.name;
	numTriangles = group.numTriangles;

	triangles.resize(group.triangles.size());
	triangles.set(group.triangles);

	rangeMaterials.resize(group.rangeMaterials.size());
	rangeMaterials.set(group.rangeMaterials);

	materials.resize(group.materials.size());
	materials.set(group.materials);
}

MeshGroup& MeshGroup::operator=(const OBJGroup &group) {
	set(group);
	return (*this);
}
