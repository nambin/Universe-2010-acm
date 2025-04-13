#include <stdafx.h>
#include <Principle/Geometry/Surface/OBJStruct.h>

using namespace geometryspace;
using namespace std;

void OBJGroup::glDraw(const OBJStruct &objstruct, int mode) const {

	const vector<Vector3<float>> &posVerts = objstruct.posVerts;

	for(unsigned int i=0;i<materials.size();++i) {

		if (mode & OBJStruct::OBJ_MATERIAL) {
			const OBJMaterial &material = objstruct.materials[materials[i]];

			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material.ambient.getptr());
			glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material.diffuse.getptr());
			glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material.specular.getptr());
			glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, material.shininess);
		}

		if (mode & OBJStruct::OBJ_COLOR) {
			const OBJMaterial &material = objstruct.materials[materials[i]];
			glColor3fv(material.diffuse.getptr());
		}
		if (mode & OBJStruct::OBJ_TEXTURE) {
			const OBJMaterial &material = objstruct.materials[materials[i]];

			if(material.map_Kd_texture > -1) {
				glEnable(GL_TEXTURE_2D);
				glBindTexture(GL_TEXTURE_2D, static_cast<GLuint>(material.map_Kd_texture));
			}
			else 
				glDisable(GL_TEXTURE_2D);
		}			

		glBegin(GL_TRIANGLES);

		for(int j=rangeMaterials[i];j<rangeMaterials[i+1];++j) {
			const OBJTriangle &triangle = objstruct.triangles[triangles[j]];

			if(mode & OBJStruct::OBJ_FLAT)
				glNormal3fv(objstruct.faceNormals[triangle.fIndex].getptr());

			for(int k=0;k<3;++k) {
				if(mode & OBJStruct::OBJ_SMOOTH)
					glNormal3fv(objstruct.normals[triangle.nIndex[k]].getptr());
				if(mode & OBJStruct::OBJ_TEXTURE)
					glTexCoord2fv(objstruct.texCoords[triangle.tIndex[k]].getptr());

				glVertex3fv(objstruct.posVerts[triangle.vIndex[k]].getptr());
			}			
		}

		glEnd();
	}
}

void OBJStruct::glDraw() const {

	// do a bit of warning 
	if (drawMode & OBJStruct::OBJ_FLAT && faceNormals.empty()) {
		std::cout << "OBJStruct::glDraw() warning _ flat render drawMode requested with no facet normals defined" << std::endl;
		drawMode &= ~OBJStruct::OBJ_FLAT;
	}
	if (drawMode & OBJStruct::OBJ_SMOOTH && normals.empty()) {
		std::cout << "OBJStruct::glDraw() warning: smooth render drawMode requested with no normals defined" << std::endl;
		drawMode &= ~OBJStruct::OBJ_SMOOTH;
	}
	if (drawMode & OBJStruct::OBJ_TEXTURE && texCoords.empty()) {
		std::cout << "OBJStruct::glDraw() warning: texture render drawMode requested with no texture coordinates defined" << std::endl;
		drawMode &= ~OBJStruct::OBJ_TEXTURE;
	}
	if (drawMode & OBJStruct::OBJ_FLAT && drawMode & OBJStruct::OBJ_SMOOTH) {
		std::cout << "OBJStruct::glDraw() warning: flat render drawMode requested and smooth render drawMode requested (using smooth)" << std::endl;
		drawMode &= ~OBJStruct::OBJ_FLAT;
	}
	if (drawMode & OBJStruct::OBJ_COLOR && materials.empty()) {
		std::cout << "OBJStruct::glDraw() warning: color render drawMode requested with no materials defined" << std::endl;
		drawMode &= ~OBJStruct::OBJ_COLOR;
	}
	if (drawMode & OBJStruct::OBJ_MATERIAL && materials.empty()) {
		std::cout << "OBJStruct::glDraw() warning: material render drawMode requested with no materials defined" << std::endl;
		drawMode &= ~OBJStruct::OBJ_MATERIAL;
	}
	if (drawMode & OBJStruct::OBJ_COLOR && drawMode & OBJStruct::OBJ_MATERIAL) {
		std::cout << "OBJStruct::glDraw() warning: color and material render drawMode requested using only material drawMode" << std::endl;
		drawMode &= ~OBJStruct::OBJ_COLOR;
	}

	if (drawMode & OBJStruct::OBJ_COLOR)			glEnable(GL_COLOR_MATERIAL);
	else if (drawMode & OBJStruct::OBJ_MATERIAL)	glDisable(GL_COLOR_MATERIAL);

	OBJGroup *group = groups;
	while(group) {
		group->glDraw((*this), drawMode);
		group = group->next;
	}

	glDisable(GL_TEXTURE_2D);
}
