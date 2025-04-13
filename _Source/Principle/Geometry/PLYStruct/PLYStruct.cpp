#include <stdafx.h>
#include <sstream>
#include <fstream>
#include <Principle/VirtualSpace/RenderingPipeline/RenderingPipelineSetting.h>
#include <Principle/Geometry/Surface/PLYStruct.h>
#include <Principle/Geometry/PLYLibrary/PLYLoader.h>

#pragma warning(disable:4800)

using namespace plyloader;

namespace geometryspace {

void PLYStruct::init() {
	numVertices = numTriangles = 0;
	vertices.clear();
	triangles.clear();
}

void PLYStruct::loadFromPLY(const char* const _filename) {

	std::ifstream is(_filename);
	if(!is.is_open()) {
		basis::throwError("PLY file open error : %s\n", _filename);
		return;
	}
	is.close();
	
	init();

	//////////////////////////////////////////////////////////////////////////
	// Load using plyLoader
	int nverts, nfaces;
	Vertex **vlist;
	Face **flist;

	int per_vertex_color = 0, has_normals = 0;

	FILE *fp = NULL;
	fopen_s(&fp, _filename, "rb");
	plyLoader(fp, nverts, nfaces, vlist, flist, per_vertex_color, has_normals);

	//////////////////////////////////////////////////////////////////////////
	// Import Vertices & Faces
	PLYVertex plyvert;
	for(int i=0;i<nverts;++i) {
		memcpy(&plyvert.pos[0], &((*vlist[i]).x), sizeof(float) * 7 + sizeof(unsigned char) * 3);				
		vertices.push_back(plyvert);
	}
	numVertices = static_cast<int>(vertices.size());

	PLYTriangle plytri;
	for(int i=0;i<nfaces;++i) {
		memcpy(plytri.vIndex, (*flist[i]).verts, sizeof(int) * 3);
		triangles.push_back(plytri);

		for(int idx=3;idx<(*flist[i]).nverts;++idx) {
			plytri.vIndex[1] = plytri.vIndex[2];
			plytri.vIndex[2] = (*flist[i]).verts[idx];
			triangles.push_back(plytri);
		}
	}
	numTriangles = static_cast<int>(triangles.size());
//	for(unsigned int i=0;i<triangles.size();++i)
//		std::swap(triangles[i].vIndex[0], triangles[i].vIndex[1]);	

	hasVertexColor = static_cast<bool>(per_vertex_color);
	hasNormals = static_cast<bool>(has_normals);

	//////////////////////////////////////////////////////////////////////////
	// Freeing
	for(int i=0;i<nverts;++i) 
		free(vlist[i]);
	free(vlist);

	for(int i=0;i<nfaces;++i) {
		free((*flist[i]).verts);
		free(flist[i]);
	}
	free(flist);

	fclose(fp);
}

void PLYStruct::glDrawSimpleMesh(const RGBA &color) const {

	glEnable(GL_COLOR_MATERIAL);
	glColor4f(color[0], color[1], color[2], color[3]);

	glBegin(GL_TRIANGLES);
	if(hasNormals) {
		for(int i=0;i<numTriangles;++i) {				
			for(int k=0;k<3;++k) {
				const PLYVertex &vert = vertices[triangles[i].vIndex[k]];
				glNormal3f(vert.normal[0], vert.normal[1], vert.normal[2]);
				glVertex3f(vert.pos[0], vert.pos[1], vert.pos[2]);
			}					
		}		
	}
	else {		
		glNormal3f(1.0f, 0.0f, 0.0f);
		for(int i=0;i<numTriangles;++i) 
			for(int k=0;k<3;++k) {
				const PLYVertex &vert = vertices[triangles[i].vIndex[k]];
				glVertex3f(vert.pos[0], vert.pos[1], vert.pos[2]);
			}
	}	
	glEnd();
}

void PLYStruct::glDrawVertexColoredMesh() const {

	glDisable(GL_COLOR_MATERIAL);
	glDisable(GL_LIGHTING);

	glBegin(GL_TRIANGLES);
	glNormal3f(1.0f, 0.0f, 0.0f);
	for(int i=0;i<numTriangles;++i) {
		for(int k=0;k<3;++k) {
			const PLYVertex &vert = vertices[triangles[i].vIndex[k]];
			glColor3ub(vert.rgb[0], vert.rgb[1], vert.rgb[2]);
			glVertex3f(vert.pos[0], vert.pos[1], vert.pos[2]);
		}
	}
	glEnd();
}

std::ostream &operator<<(std::ostream &os, const PLYVertex &v) {
	os << "XYZ : " << v.pos << std::endl;
	os << "RGB : " << v.rgb << std::endl;
	os << "NORMAL : " << v.normal << std::endl;
	return os;
}

std::ostream &operator<<(std::ostream &os, const PLYTriangle &f) {	
	os << "VertIdx : ";
	for(int i=0;i<3;++i) 
		os << f.vIndex[i] << " ";
	os << std::endl;
	return os;
}

};			// geometryspace