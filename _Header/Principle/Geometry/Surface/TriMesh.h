#ifndef _GEOMETRYSPACE_TRIMESH_H_
#define _GEOMETRYSPACE_TRIMESH_H_

#include <sstream>
#include <fstream>
#include "../../../Basis/StringLibrary/StringLibrary.h"
#include "../../Arithmetic/Vector/Vector2.h"
#include "../../Arithmetic/Vector/Vector3.h"
#include "../../VirtualSpace/Material/Material.h"
#include "TriMeshConnectivity.h"
#include "OBJStruct.h"
#include "PLYStruct.h"
#include "TriMesh_Simple.h"

namespace geometryspace {

using namespace std;
using basis::RGBA;
using arithmetic::Vector2;
using arithmetic::Vector3;
using virtualspace::Material;

struct PRINCIPLEDLL MeshMaterial {
public : 
	MeshMaterial() { init(); }
	
	void init();	
	void set(const MeshMaterial &mat);
	void set(const OBJMaterial &mat);
	void set_DefaultMaterial();

	MeshMaterial& operator=(const OBJMaterial &mat);

public :
	Material	material;
	string		map_Kd_filename;
	int			map_Kd_texture;
};

template<typename POSVECTOR> class TriMesh;
struct PRINCIPLEDLL MeshGroup {
public : 
	MeshGroup(std::string _name = std::string("")) { init(_name); }
	
	void init(std::string _name);
	void set(const MeshGroup &group);
	void set(const OBJGroup &group);

	MeshGroup& operator=(const OBJGroup &group);

	template<typename POSVECTOR>
	void glDraw(const TriMesh<POSVECTOR> &mesh, int mode) const;

public : 
	std::string		name;          
	int				numTriangles;
	Array<int>		triangles;

	Array<int>		rangeMaterials;
	Array<int>		materials;	
};

template<typename POSVECTOR>
class TriMesh : public TriMeshConnectivity {
public : 
	static const int MESH_NONE		= 0;			// render with only vertices 
	static const int MESH_FLAT		= (1 << 0);		// render with facet normals 
	static const int MESH_SMOOTH	= (1 << 1);		// render with vertex normals 
	static const int MESH_TEXTURE	= (1 << 2);     // render with texture coords 
	static const int MESH_COLOR		= (1 << 3);		// render with colors 
	static const int MESH_MATERIAL	= (1 << 4);		// render with materials 

public:
	typedef typename POSVECTOR::value_type real_type;
	typedef typename POSVECTOR::value_type T;

	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	TriMesh();	
	template<typename _POSVECTOR> TriMesh(const TriMesh<_POSVECTOR> &mesh);

	//////////////////////////////////////////////////////////////////////////
	// Initialization
private :
	void init();

public : 
	template<typename _POSVECTOR> void set(const TriMesh<_POSVECTOR> &mesh);
	void set(const OBJStruct &obj);
	void set(const PLYStruct &ply);
	template<typename CONTAINER> void set(const TriMesh_Simple<POSVECTOR, CONTAINER>&);
	void release();

	//////////////////////////////////////////////////////////////////////////
	//	Loader
	void load(const char * const filename);
	void loadFromOBJ(const char * const filename);
	void loadFromPLY(const char * const filename);

	void makeDefaultGroup();

	//////////////////////////////////////////////////////////////////////////
	//	glDraw
	void glDraw() const;
	void glDraw(int mode) const;

	void glDrawLineStrip(float width, const RGBA &color = RGBA::white())		const;
	void glDrawLineStrip2D(float width, const RGBA &color = RGBA::white())		const;
	void glDrawSmoothSurface(const Material& mat = RGBA::white())				const;
	void glDrawSimpleMesh(const RGBA &color = RGBA::white())					const;	
	void glDrawTriangle(int idx, const RGBA &color = RGBA::red())				const;	
	void glDrawTriangle2D(int idx, const RGBA &color = RGBA::red())				const;	
	void glDrawEdge2D(int idx, float width, const RGBA& color = RGBA::red())	const;
	void glDrawVertex(int idx, float size, const RGBA& color = RGBA::red())		const;
	void glDrawVertex2D(int idx, float size, const RGBA& color = RGBA::red())	const;
	void glDrawNormalMap()														const;
	void glDrawFaceNormalMap()													const;

	void glDrawBoundaryVertex2D(float size, const RGBA& c0 = RGBA::red()) const;
	void glDrawV2V_2D(int iv, int nth, float size, const RGBA& c0 = RGBA::red(), const RGBA& c1 = RGBA::blue()) const;
	void glDrawE2V_2D(int iv, int nth, float size, const RGBA& c0 = RGBA::red(), const RGBA& c1 = RGBA::blue()) const;
	void glDrawF2V_2D(int iv, int nth, float size, const RGBA& c0 = RGBA::red(), const RGBA& c1 = RGBA::blue()) const;

	//////////////////////////////////////////////////////////////////////////
	// Manipulation
	void setDrawMode(int mode) const { drawMode = mode; }
	void makeFaceNormals();
	void makeVertexNormals();

	void changeTriangleOrder();

	void unitize();
	void translate(const POSVECTOR &trans);
	void scale(real_type s);	

public :
	//////////////////////////////////////////////////////////////////////////
	// Numbers
	int							numNormals;
	int							numTexCoords;
	int							numMaterials;
	int							numGroups;

	//////////////////////////////////////////////////////////////////////////
	// Geometry	
	Array2<int>					n2f;				// [numFaces*3]		normals indices of each face
	Array2<int>					t2f;				// [numFaces*3]		texCoords indices of each face
	
	Array<POSVECTOR>			posVerts;			// [numVertices]
	Array<Vector3<real_type>>	normals;			// [numNormals]
	Array<Vector3<real_type>>	faceNormals;		// [numFaces]
	Array<Vector2<real_type>>	texCoords;			// [numTexCoords]
	
	//////////////////////////////////////////////////////////////////////////
	// Material
	Array<MeshMaterial>			materials;
	Array<MeshGroup>			groups;	

	//////////////////////////////////////////////////////////////////////////
	// Rendering
	mutable int					drawMode;
};

template<typename POSVECTOR>
TriMesh<POSVECTOR>::TriMesh() : TriMeshConnectivity() {
	init();
}

template<typename POSVECTOR> template<typename _POSVECTOR> 
TriMesh<POSVECTOR>::TriMesh(const TriMesh<_POSVECTOR> &mesh) : TriMeshConnectivity() {
	init(); 
	set(mesh);
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::init() {
	numNormals = numTexCoords = numMaterials = numGroups = 0;

	n2f.clear(); t2f.clear();	
	posVerts.clear(); normals.clear(); faceNormals.clear(); texCoords.clear();
	materials.clear(); groups.clear();

	drawMode = (MESH_FLAT | MESH_SMOOTH | MESH_TEXTURE | MESH_COLOR | MESH_MATERIAL);
}

template<typename POSVECTOR> template<typename _POSVECTOR>
void TriMesh<POSVECTOR>::set(const TriMesh<_POSVECTOR> &mesh) {

	TriMeshConnectivity::release();
	release();

	TriMeshConnectivity::set(static_cast<const TriMeshConnectivity&>(mesh));

	//////////////////////////////////////////////////////////////////////////
	// Numbers
	numNormals = mesh.numNormals;
	numTexCoords = mesh.numTexCoords;
	numMaterials = mesh.numMaterials;
	numGroups = mesh.numGroups;

	//////////////////////////////////////////////////////////////////////////
	// Geometry
	n2f.resize_memcpy(mesh.n2f);
	t2f.resize_memcpy(mesh.t2f);

	posVerts.resize(mesh.posVerts.size());
	normals.resize(mesh.normals.size());
	faceNormals.resize(mesh.faceNormals.size());
	texCoords.resize(mesh.texCoords.size());
	
	for(size_t i=0;i<posVerts.size();++i)		posVerts[i] = mesh.posVerts[i];
	for(size_t i=0;i<normals.size();++i)		normals[i] = mesh.normals[i];
	for(size_t i=0;i<faceNormals.size();++i)	faceNormals[i] = mesh.faceNormals[i];
	for(size_t i=0;i<texCoords.size();++i)		texCoords[i] = mesh.texCoords[i];

	//////////////////////////////////////////////////////////////////////////
	// Material
	materials.resize_set(mesh.materials);
	groups.resize_set(mesh.groups);
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::set(const OBJStruct &obj) {
	TriMeshConnectivity::release();
	release();

	numVertices = obj.numVertices;
	numFaces = obj.numTriangles;
	numNormals = obj.numNormals;
	numTexCoords = obj.numTexCoords;
	numMaterials = obj.numMaterials;
	numGroups = obj.numGroups;

	v2f.resize(numFaces, 3);
	for(int i=0;i<numFaces;++i) for(int j=0;j<3;++j) v2f(i,j) = obj.triangles[i].vIndex[j];

	n2f.resize(numFaces, 3);
	for(int i=0;i<numFaces;++i) for(int j=0;j<3;++j) n2f(i,j) = obj.triangles[i].nIndex[j];

	t2f.resize(numFaces, 3);
	for(int i=0;i<numFaces;++i) for(int j=0;j<3;++j) t2f(i,j) = obj.triangles[i].tIndex[j];

	posVerts.resize(numVertices);
	for(int i=0;i<numVertices;++i) posVerts[i] = obj.posVerts[i];

	normals.resize(numNormals);
	for(int i=0;i<numNormals;++i) normals[i] = obj.normals[i];

	texCoords.resize(numTexCoords);
	for(int i=0;i<numTexCoords;++i) texCoords[i] = obj.texCoords[i];

	materials.resize(numMaterials);
	for(int i=0;i<numMaterials;++i) materials[i] = obj.materials[i];

	groups.resize(numGroups);
	OBJGroup *objgroup = obj.groups;
	for(int i=0;i<numGroups;++i) {		
		groups[i] = (*objgroup);
		objgroup = objgroup->next;
	}
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::set(const PLYStruct &ply) {
	TriMeshConnectivity::release();
	release();
	
	numVertices = ply.numVertices;
	numFaces = ply.numTriangles;

	v2f.resize(numFaces, 3);
	for(int i=0;i<numFaces;++i) for(int j=0;j<3;++j) v2f(i,j) = ply.triangles[i].vIndex[j];

	posVerts.resize(numVertices);
	for(int i=0;i<numVertices;++i) posVerts[i] = ply.vertices[i].pos;

	if(ply.hasNormals) {
		numNormals = numVertices;
		normals.resize(numNormals);
		for(int i=0;i<numNormals;++i) normals[i] = ply.vertices[i].normal;
	}
}

template<typename POSVECTOR> template<typename CONTAINER> 
void TriMesh<POSVECTOR>::set(const TriMesh_Simple<POSVECTOR, CONTAINER>& mesh) {
	numVertices = mesh.numVertices;
	numFaces = mesh.numFaces;

	posVerts.resize(numVertices);
	v2f.resize(numFaces, 3);

	for(int i=0;i<numFaces;++i) for(int j=0;j<3;++j)	v2f(i,j) = mesh.v2f[i*3+j];
	for(int i=0;i<numVertices;++i)						posVerts[i] = mesh.posVerts[i];
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::release() {
	init();
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::load(const char * const filename) {

	std::string ext = basis::DataPath::getFilenameExtension(filename);

	if(ext.compare("OBJ")==0 || ext.compare("obj")==0)		return loadFromOBJ(filename);
	else if(ext.compare("PLY")==0 || ext.compare("ply")==0)	return loadFromPLY(filename);
	else
		basis::throwError("Unsupported FileFormat _ TriMesh::load(..)\n");

	return;
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::loadFromOBJ(const char* const filename) {
	OBJStruct obj;
	obj.loadFromOBJ(filename);
	set(obj);
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::loadFromPLY(const char* const filename) {
	PLYStruct ply;
	ply.loadFromPLY(filename);
	set(ply);
	
	makeDefaultGroup();

	std::cout << "Load From PLY done " << numVertices << " " << numFaces << endl;
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::makeDefaultGroup() {
	numMaterials = 1;
	numGroups = 1;

	materials.resize(numMaterials);
	materials[0].set_DefaultMaterial();

	groups.resize(numGroups);
	groups[0].numTriangles = numFaces;
	groups[0].triangles.resize(numFaces);
	for(size_t i=0;i<groups[0].triangles.size();++i)
		groups[0].triangles[i] = static_cast<int>(i);
	groups[0].materials.resize(1);
	groups[0].rangeMaterials.resize(2);
	groups[0].materials[0] = 0;
	groups[0].rangeMaterials[0] = 0;
	groups[0].rangeMaterials[1] = numFaces;
}

template<typename POSVECTOR>
void MeshGroup::glDraw(const TriMesh<POSVECTOR> &mesh, int mode) const {

	const Array<POSVECTOR> &posVerts = mesh.posVerts;

	for(unsigned int i=0;i<materials.size();++i) {

		if (mode & TriMesh<POSVECTOR>::MESH_MATERIAL) {
			const Material &material = mesh.materials[materials[i]].material;

			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material.getAmbient().getptr());
			glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material.getDiffuse().getptr());
			glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material.getSpecular().getptr());
			glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, material.getShininess());
		}

		if (mode & TriMesh<POSVECTOR>::MESH_COLOR) {
			const Material &material = mesh.materials[materials[i]].material;
			glColor3fv(material.getDiffuse().getptr());
		}
		if (mode & TriMesh<POSVECTOR>::MESH_TEXTURE) {
			const MeshMaterial &material = mesh.materials[materials[i]];

			if(material.map_Kd_texture > -1) {
				glEnable(GL_TEXTURE_2D);
				glBindTexture(GL_TEXTURE_2D, static_cast<GLuint>(material.map_Kd_texture));
			}
			else 
				glDisable(GL_TEXTURE_2D);
		}			

		glBegin(GL_TRIANGLES);

		for(int j=rangeMaterials[i];j<rangeMaterials[i+1];++j) {
			if(mode & TriMesh<POSVECTOR>::MESH_FLAT)
				glNormal3fv(mesh.faceNormals[triangles[j]].getptr());

			for(int k=0;k<3;++k) {
				if(mode & TriMesh<POSVECTOR>::MESH_SMOOTH)
					glNormal3fv(mesh.normals[mesh.n2f(triangles[j], k)].getptr());
				if(mode & TriMesh<POSVECTOR>::MESH_TEXTURE)
					glTexCoord2fv(mesh.texCoords[mesh.t2f(triangles[j], k)].getptr());

				glVertex3fv(mesh.posVerts[mesh.v2f(triangles[j], k)].getptr());
			}			
		}

		glEnd();
	}
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::glDraw() const {
	glDraw(drawMode);
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::glDraw(int mode) const {
	
	// do a bit of warning 
	if (mode & TriMesh::MESH_FLAT && faceNormals.empty()) {
		//std::cout << "TriMesh::glDraw() warning _ flat render mode requested with no facet normals defined" << std::endl;
		mode &= ~TriMesh::MESH_FLAT;
	}
	if (mode & TriMesh::MESH_SMOOTH && normals.empty()) {
		//std::cout << "TriMesh::glDraw() warning: smooth render mode requested with no normals defined" << std::endl;
		mode &= ~TriMesh::MESH_SMOOTH;
	}
	if (mode & TriMesh::MESH_TEXTURE && texCoords.empty()) {
		//std::cout << "TriMesh::glDraw() warning: texture render mode requested with no texture coordinates defined" << std::endl;
		mode &= ~TriMesh::MESH_TEXTURE;
	}
	if (mode & TriMesh::MESH_FLAT && mode & TriMesh::MESH_SMOOTH) {
		//std::cout << "TriMesh::glDraw() warning: flat render mode requested and smooth render mode requested (using smooth)" << std::endl;
		mode &= ~TriMesh::MESH_FLAT;
	}
	if (mode & TriMesh::MESH_COLOR && materials.empty()) {
		//std::cout << "TriMesh::glDraw() warning: color render mode requested with no materials defined" << std::endl;
		mode &= ~TriMesh::MESH_COLOR;
	}
	if (mode & TriMesh::MESH_MATERIAL && materials.empty()) {
		//std::cout << "TriMesh::glDraw() warning: material render mode requested with no materials defined" << std::endl;
		mode &= ~TriMesh::MESH_MATERIAL;
	}
	if (mode & TriMesh::MESH_COLOR && mode & TriMesh::MESH_MATERIAL) {
		//std::cout << "TriMesh::glDraw() warning: color and material render mode requested using only material mode" << std::endl;
		mode &= ~TriMesh::MESH_COLOR;
	}

	if (mode & TriMesh::MESH_COLOR)			glEnable(GL_COLOR_MATERIAL);
	else if (mode & TriMesh::MESH_MATERIAL)	glDisable(GL_COLOR_MATERIAL);

	for(int i=0;i<numGroups;++i) 
		groups[i].glDraw((*this), mode);

	if(numGroups == 0)
		glDrawSimpleMesh(RGBA::white());

	glDisable(GL_TEXTURE_2D);
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::glDrawLineStrip(float width, const RGBA &color) const {

	glLineWidth(width);
	glEnable(GL_COLOR_MATERIAL);

	glColor4f(color[0], color[1], color[2], color[3]);
	glBegin(GL_LINES);
	for(int i=0;i<numEdges;++i) {		
		glVertex3f(posVerts[v2e[i*2+0]][0], posVerts[v2e[i*2+0]][1], posVerts[v2e[i*2+0]][2]);
		glVertex3f(posVerts[v2e[i*2+1]][0], posVerts[v2e[i*2+1]][1], posVerts[v2e[i*2+1]][2]);
	}
	glEnd();

	cout << numEdges << endl;
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::glDrawLineStrip2D(float width, const RGBA &color) const {

	glLineWidth(width);
	glEnable(GL_COLOR_MATERIAL);

	glColor4f(color[0], color[1], color[2], color[3]);
	glBegin(GL_LINES);
	for(int i=0;i<numEdges;++i) {		
		glVertex3f(posVerts[v2e[i*2+0]][0], posVerts[v2e[i*2+0]][1], 0);
		glVertex3f(posVerts[v2e[i*2+1]][0], posVerts[v2e[i*2+1]][1], 0);
	}
	glEnd();
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::glDrawSmoothSurface(const Material& mat = RGBA::white()) const {
	if(normals.empty()) {
		glDrawSimpleMesh(mat.getDiffuse());
		return;
	}

	glEnable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);

	glDisable(GL_COLOR_MATERIAL);
	mat.glBindMaterial();	

	glBegin( GL_TRIANGLES );

	for(int i=0;i<numFaces;++i)	{
		unsigned int i0 = v2f[i*3  ];
		unsigned int i1 = v2f[i*3+1];
		unsigned int i2 = v2f[i*3+2];

		const Vector3<real_type> *V0, *V1, *V2, *N0, *N1, *N2;
		V0 = &posVerts[v2f[i*3]  ];
		V1 = &posVerts[v2f[i*3+1]];
		V2 = &posVerts[v2f[i*3+2]];
		N0 = &normals[n2f[i*3]  ];
		N1 = &normals[n2f[i*3+1]];
		N2 = &normals[n2f[i*3+2]];

		glNormal3f((*N0)[0], (*N0)[1], (*N0)[2]);
		glVertex3f((*V0)[0], (*V0)[1], (*V0)[2]);
		glNormal3f((*N1)[0], (*N1)[1], (*N1)[2]);
		glVertex3f((*V1)[0], (*V1)[1], (*V1)[2]);
		glNormal3f((*N2)[0], (*N2)[1], (*N2)[2]);
		glVertex3f((*V2)[0], (*V2)[1], (*V2)[2]);
	}

	glEnd();
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::glDrawSimpleMesh(const RGBA &color) const {	

	glEnable(GL_COLOR_MATERIAL);
	
	glColor4f(color[0], color[1], color[2], color[3]);
	glBegin(GL_TRIANGLES);
	if(normals.size() != 0) {
		for(int i=0;i<numFaces;++i) {				
			for(int k=0;k<3;++k) {
				glNormal3f(normals[n2f[i*3+k]][0], normals[n2f[i*3+k]][1], normals[n2f[i*3+k]][2]);
				glVertex3f(posVerts[v2f[i*3+k]][0], posVerts[v2f[i*3+k]][1], posVerts[v2f[i*3+k]][2]);			
			}					
		}		
	}
	else {		
		glNormal3f(1.0f, 0.0f, 0.0f);
		for(int i=0;i<numFaces;++i) 
			for(int k=0;k<3;++k) 				
				glVertex3f(posVerts[v2f[i*3+k]][0], posVerts[v2f[i*3+k]][1], posVerts[v2f[i*3+k]][2]);			
	}	
	glEnd();
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::glDrawTriangle(int idx, const RGBA &color) const {
	
	glEnable(GL_COLOR_MATERIAL);

	glColor4f(color[0], color[1], color[2], color[3]);
	glBegin(GL_TRIANGLES);
	if(normals.size() != 0) {
		for(int k=0;k<3;++k) {
			glNormal3f(normals[n2f[idx*3+k]][0], normals[n2f[idx*3+k]][1], normals[n2f[idx*3+k]][2]);
			glVertex3f(posVerts[v2f[idx*3+k]][0], posVerts[v2f[idx*3+k]][1], posVerts[v2f[idx*3+k]][2]);			
		}					
	}
	else {
		glNormal3f(1.0f, 0.0f, 0.0f);
		for(int k=0;k<3;++k) 
			glVertex3f(posVerts[v2f[idx*3+k]][0], posVerts[v2f[idx*3+k]][1], posVerts[v2f[idx*3+k]][2]);			
	}
	glEnd();
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::glDrawTriangle2D(int idx, const RGBA &color) const {

	glEnable(GL_COLOR_MATERIAL);

	glColor4f(color[0], color[1], color[2], color[3]);
	glBegin(GL_TRIANGLES);
		for(int k=0;k<3;++k) 
			glVertex3f(posVerts[v2f(idx,k)][0], posVerts[v2f(idx,k)][1], 0.001f);
	glEnd();
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::glDrawEdge2D(int idx, float size, const RGBA& color) const {
	glLineWidth(size);
	glEnable(GL_COLOR_MATERIAL);
	glColor4f(color[0], color[1], color[2], color[3]);
	glBegin(GL_LINES);
		glVertex3f(posVerts[v2e(idx,0)][0], posVerts[v2e(idx,0)][1], 0.001f);
		glVertex3f(posVerts[v2e(idx,1)][0], posVerts[v2e(idx,1)][1], 0.001f);
	glEnd();
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::glDrawVertex(int idx, float size, const RGBA& color) const {
	glPointSize(size);
	glEnable(GL_COLOR_MATERIAL);
	glColor4f(color[0], color[1], color[2], color[3]);
	glBegin(GL_POINTS);
	glVertex3f(posVerts[idx][0], posVerts[idx][1], posVerts[idx][2]);
	glEnd();
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::glDrawVertex2D(int idx, float size, const RGBA& color) const {
	glPointSize(size);
	glEnable(GL_COLOR_MATERIAL);
	glColor4f(color[0], color[1], color[2], color[3]);
	glBegin(GL_POINTS);
		glVertex3f(posVerts[idx][0], posVerts[idx][1], 0.001f);
	glEnd();
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::glDrawNormalMap() const {

	glEnable(GL_COLOR_MATERIAL);
	
	if(normals.size() == 0 && faceNormals.size() == 0) {
		glDrawSimpleMesh();
		return;
	}
	else if(normals.size() == 0 && faceNormals.size() != 0) {
		glDrawFaceNormalMap();
		return;
	}

	glBegin(GL_TRIANGLES);
	for(int i=0;i<numFaces;++i) {				
		for(int k=0;k<3;++k) {
			glColor4f(normals[n2f[i*3+k]][0], normals[n2f[i*3+k]][1], normals[n2f[i*3+k]][2], 1.0f);
			glNormal3f(normals[n2f[i*3+k]][0], normals[n2f[i*3+k]][1], normals[n2f[i*3+k]][2]);
			glVertex3f(posVerts[v2f[i*3+k]][0], posVerts[v2f[i*3+k]][1], posVerts[v2f[i*3+k]][2]);			
		}					
	}		
	glEnd();
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::glDrawFaceNormalMap() const {

	glEnable(GL_COLOR_MATERIAL);

	if(faceNormals.size() == 0) {
		glDrawSimpleMesh();
		return;
	}

	glBegin(GL_TRIANGLES);
	for(int i=0;i<numFaces;++i) {				
		glColor4f(faceNormals[i][0], faceNormals[i][1], faceNormals[i][2], 1.0f);
		glNormal3f(faceNormals[i][0], faceNormals[i][1], faceNormals[i][2]);
		for(int k=0;k<3;++k) 		
			glVertex3f(posVerts[v2f[i*3+k]][0], posVerts[v2f[i*3+k]][1], posVerts[v2f[i*3+k]][2]);			
	}		
	glEnd();
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::glDrawBoundaryVertex2D(float size, const RGBA& c0) const {
	for(int i=0;i<numVertices;++i) { if(isBdryVertex[i]) { glDrawVertex2D(i, size, c0); }}
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::glDrawV2V_2D(int iv, int nth, float size, const RGBA& c0, const RGBA& c1) const {
	glDrawVertex2D(iv, size, c0);
	if(nth < (int)v2v[iv].size()) glDrawVertex2D(v2v[iv][nth], size, c1);
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::glDrawE2V_2D(int iv, int nth, float size, const RGBA& c0, const RGBA& c1) const {
	glDrawVertex2D(iv, size, c0);
	if(nth < (int)e2v[iv].size()) glDrawEdge2D(e2v[iv][nth], size, c1);
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::glDrawF2V_2D(int iv, int nth, float size, const RGBA& c0, const RGBA& c1) const {
	glDrawVertex2D(iv, size, c0);
	if(nth < (int)f2v[iv].size()) glDrawTriangle2D(f2v[iv][nth], c1);
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::makeFaceNormals() { 

	if(static_cast<int>(faceNormals.size()) == numFaces)
		return;

	faceNormals.resize(numFaces);

	Vector3<double> n, v0, v1, v2;
	for(int i=0;i<numFaces;++i) {		
		v0 = posVerts[v2f(i,0)];
		v1 = posVerts[v2f(i,1)];
		v2 = posVerts[v2f(i,2)];

		n = arithmetic::unitVector(arithmetic::cross(v1-v0, v2-v0));
		faceNormals[i] = n;
	}
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::makeVertexNormals() {
	normals.resize(posVerts.size());
	numNormals = normals.size();

	if(faceNormals.size() != numFaces)	
		makeFaceNormals();

	Array<vector<unsigned int>> f2v(numVertices);
	for(int i=0;i<numFaces;i++)	{
		f2v[v2f[i*3  ]].push_back(i);
		f2v[v2f[i*3+1]].push_back(i);
		f2v[v2f[i*3+2]].push_back(i);
	}
	
	for(size_t i=0;i<normals.size();++i) {
		normals[i] = Vector3<T>(0.0,0.0,0.0);
		
		for(size_t k=0;k<f2v[i].size();++k) 
			normals[i] += faceNormals[f2v[i][k]];

		normals[i].normalize();
	}

	n2f.resize_memcpy(v2f);
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::changeTriangleOrder() {
	for(size_t i=0;i<v2f.dimX();++i) 
		std::swap(v2f(i,1), v2f(i,2));
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::unitize() {

	real_type maxx, minx, maxy, miny, maxz, minz;

	// get the max/min
	maxx = minx = posVerts[0][0];
	maxy = miny = posVerts[0][1];
	maxz = minz = posVerts[0][2];

	for(int i = 1; i < numVertices; ++i) {
		if (maxx < posVerts[i][0])	maxx = posVerts[i][0];
		if (minx > posVerts[i][0])	minx = posVerts[i][0];

		if (maxy < posVerts[i][1])	maxy = posVerts[i][1];
		if (miny > posVerts[i][1])	miny = posVerts[i][1];

		if (maxz < posVerts[i][2])	maxz = posVerts[i][2];
		if (minz > posVerts[i][2])	minz = posVerts[i][2];
	}

	// calculate model width, height, and depth 
	real_type w = abs(maxx-minx);
	real_type h = abs(maxy-miny);
	real_type d = abs(maxz-minz);

	// calculate center of the model
	real_type cx = (maxx + minx) * 0.5f;
	real_type cy = (maxy + miny) * 0.5f;
	real_type cz = (maxz + minz) * 0.5f;

	// calculate unitizing scale factor 
	real_type s = 1.0f / basis::MAX3(w, h, d);

	// translate around center then scale 
	translate(-POSVECTOR(cx, cy, cz));
	scale(s);

//	translate(-POSVECTOR(minx, miny, minz));
//	scale(s);	
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::translate(const POSVECTOR &trans) {
	for(int i = 0; i < numVertices; ++i) 
		posVerts[i] += trans;
}

template<typename POSVECTOR>
void TriMesh<POSVECTOR>::scale(real_type s)	{
	for(int i = 0; i < numVertices; ++i) {
		posVerts[i][0] *= s;
		posVerts[i][1] *= s;
		posVerts[i][2] *= s;
	}
}


};		// namespace geometryspace

#endif