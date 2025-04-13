#include <stdafx.h>
#include <sstream>
#include <fstream>
#include <Principle/Geometry/Surface/OBJStruct.h>

using namespace geometryspace;

void OBJStruct::loadFromOBJ(const char* const _filename) {

	std::ifstream is(_filename);
	if(!is.is_open()) {
		basis::throwError("OBJ file open error : %s\n", _filename);
		return;
	}
	is.close();

	init();
	filename = string(_filename);

	// Set Default Material
	OBJMaterial material;
	material.set_DefaultMaterial();
	materials.push_back(material);

	loadPass(_filename);
	setFaceNormals();

	std::cout << "Load OBJ File : " << filename.c_str() << std::endl;
	std::cout << "Number of Vertices : " << numVertices << std::endl;
	std::cout << "Number of Faces : " << numTriangles << std::endl;
}

void OBJStruct::readMTL(const char* const name) {

	std::string dir = basis::DataPath::getDirectory(filename);
	std::string mtlfilename = dir + name;

	std::ifstream is(mtlfilename.c_str());
	if(!is.is_open()) {
		cout << "MTL file open error : " << mtlfilename.c_str() << endl;
		return;
	}

	float d;
	OBJMaterial material;

	int idxMaterial = 0;

	std::istringstream iss;
	std::string str;
	while(getline(is, str)) {
		iss.clear(); iss.str(str);
		iss >> str;

		if(str.compare("#")==0)
			continue;

		else if(str.compare("newmtl")==0) {
			++idxMaterial;
			material.init();

			iss >> str;			
			material.name = str;
			materials.push_back(material);
		}
		else if(str.compare("Ka")==0) 
			iss >> materials[idxMaterial].ambient[0] >> materials[idxMaterial].ambient[1] >> materials[idxMaterial].ambient[2];
		else if(str.compare("Kd")==0) 
			iss >> materials[idxMaterial].diffuse[0] >> materials[idxMaterial].diffuse[1] >> materials[idxMaterial].diffuse[2];
		else if(str.compare("Ks")==0)
			iss >> materials[idxMaterial].specular[0] >> materials[idxMaterial].specular[1] >> materials[idxMaterial].specular[2];
		else if(str.compare("d")==0 || str.compare("Tr")==0) {
			iss >> d;
			materials[idxMaterial].ambient[3] = materials[idxMaterial].diffuse[3] = materials[idxMaterial].specular[3] = d;
		}
		else if(str.compare("Ns")==0) {
			iss >> materials[idxMaterial].shininess;
			materials[idxMaterial].shininess *= (128.0f / 1000.0f);
		}
		else if(str.compare("illum")==0) 
			iss >> materials[idxMaterial].illum;
		else if(str.compare("map_Kd")==0) {
			iss >> materials[idxMaterial].map_Kd;
			materials[idxMaterial].map_Kd_texture = static_cast<int>(makeTexture(materials[idxMaterial].map_Kd.c_str()));
		}
	}

	numMaterials = static_cast<int>(materials.size());
//	std::cout << "Number of Materials _ " << numMaterials << std::endl;
//	for(int i=0;i<numMaterials;++i) 
//		std::cout << materials[i] << std::endl;
}

const unsigned int OBJStruct::makeTexture(const char* const name) {

	std::string dir = basis::DataPath::getDirectory(filename);
	std::string imagefilename = dir + name;
	
	basis::Image textureImage;
	textureImage.load(imagefilename.c_str());

	GLuint texture;
	glEnable(GL_TEXTURE_2D);
	glGenTextures((GLsizei)1, &texture);
	glBindTexture(GL_TEXTURE_2D, texture);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

	glTexImage2D(GL_TEXTURE_2D, 0, 3, textureImage.getWidth(), textureImage.getHeight(), 0, GL_RGB, GL_FLOAT, textureImage.getPixels());
	
	return texture;
}

void OBJStruct::loadPass(const char* const filename) {

	int numvertices=0, numnormals=0, numtexcoords=0;

	int material = 0;			// current material 
	OBJGroup* group;			// current group pointer 
	OBJTriangle triangle;
	Vector3<float> vvs;
	Vector3<float> nns;
	Vector2<float> tts;

	std::ifstream is(filename);
	
	// make a default group 
	group = addGroup("default");
	
	std::istringstream iss, iss_;
	std::string str;

	while(std::getline(is, str) ) {
		iss.clear(); iss.str(str);
		iss >> str;

		if(str.compare("#")==0) 
			continue;

		else if(str.compare("mtllib")==0) {
			iss >> str;
			mtllibname = str;
			readMTL(mtllibname.c_str());
		}

		else if(str.compare("v")==0) {			
			iss >> vvs[0] >> vvs[1] >> vvs[2];
			posVerts.push_back(vvs);
			++numvertices;
		}
		else if(str.compare("vn")==0) {
			iss >> nns[0] >> nns[1] >> nns[2] ;
			normals.push_back(nns);
			++numnormals;
		}
		else if(str.compare("vt")==0) {			
			iss >> tts[0] >> tts[1] ;
			texCoords.push_back(tts);
			++numtexcoords;
		}

		else if(str.compare("usemtl")==0) {
			iss >> str;
			material = findMaterial(str.c_str());
			group->setMaterial(material);			
		}
		else if(str.compare("g")==0) {
			iss >> str;
			group = addGroup(str.c_str());
			group->setMaterial(material);			
		}

		else if(str.compare("f")==0) {						

			int vvi, tti, nni;
			unsigned int _count = 0;

			iss >> str;
			if(str.find("//") != std::string::npos) {
				for(int i=0;i<3;++i) {
					if(i != 0)
						iss >> str;

					basis::StringLibrary::replace(str, "//", " ");
					iss_.clear(); iss_.str(str);
					iss_ >> vvi >> nni;

					triangle.vIndex[i] = vvi < 0 ? vvi + numvertices : vvi-1;
					triangle.nIndex[i] = nni < 0 ? nni + numnormals : nni-1;					
				}

				triangle.fIndex = static_cast<int>(triangles.size());
				triangles.push_back(triangle);		
				group->pushTriangle(triangle.fIndex);	

				while(iss >> str) {
					basis::StringLibrary::replace(str, "//", " ");
					iss_.clear(); iss_.str(str);
					if(!(iss_ >> vvi >> nni))
						break;
					
					triangle.vIndex[0] = triangles[triangle.fIndex].vIndex[0];
					triangle.nIndex[0] = triangles[triangle.fIndex].nIndex[0];
					triangle.vIndex[1] = triangles[triangle.fIndex].vIndex[2];
					triangle.nIndex[1] = triangles[triangle.fIndex].nIndex[2];
					triangle.vIndex[2] = vvi < 0 ? vvi + numvertices : vvi-1;
					triangle.nIndex[2] = nni < 0 ? nni + numnormals : nni-1;	

					++triangle.fIndex;
					triangles.push_back(triangle);		
					group->pushTriangle(triangle.fIndex);	
				}
			}
			else if(_count = basis::StringLibrary::count(str, "/") == 2) {
				for(int i=0;i<3;++i) {	
					if(i != 0)
						iss >> str;

					basis::StringLibrary::replace(str, "/", " ");
					iss_.clear(); iss_.str(str);
					iss_ >> vvi >> tti >> nni;
					triangle.vIndex[i] = vvi < 0 ? vvi + numvertices : vvi-1;
					triangle.tIndex[i] = tti < 0 ? tti + numtexcoords : tti-1;
					triangle.nIndex[i] = nni < 0 ? nni + numnormals : nni-1;
				}

				triangle.fIndex = static_cast<int>(triangles.size());
				triangles.push_back(triangle);		
				group->pushTriangle(triangle.fIndex);	

				while(iss >> str) {
					basis::StringLibrary::replace(str, "/", " ");
					iss_.clear(); iss_.str(str);
					if(!(iss_ >> vvi >> tti >> nni))
						break;

					triangle.vIndex[0] = triangles[triangle.fIndex].vIndex[0];
					triangle.nIndex[0] = triangles[triangle.fIndex].nIndex[0];
					triangle.tIndex[0] = triangles[triangle.fIndex].tIndex[0];
					triangle.vIndex[1] = triangles[triangle.fIndex].vIndex[2];
					triangle.nIndex[1] = triangles[triangle.fIndex].nIndex[2];
					triangle.tIndex[1] = triangles[triangle.fIndex].tIndex[2];
					triangle.vIndex[2] = vvi < 0 ? vvi + numvertices : vvi-1;
					triangle.nIndex[2] = nni < 0 ? nni + numnormals : nni-1;	
					triangle.tIndex[2] = tti < 0 ? tti + numtexcoords : tti-1;

					++triangle.fIndex;
					triangles.push_back(triangle);		
					group->pushTriangle(triangle.fIndex);	
				}
			}
			else if(_count == 1) {
				for(int i=0;i<3;++i) {	
					if(i != 0)
						iss >> str;

					basis::StringLibrary::replace(str, "/", " ");
					iss_.clear(); iss_.str(str);
					iss_ >> vvi >> tti;
					triangle.vIndex[i] = vvi < 0 ? vvi + numvertices : vvi-1;
					triangle.tIndex[i] = tti < 0 ? tti + numtexcoords : tti-1;					
				}

				triangle.fIndex = static_cast<int>(triangles.size());
				triangles.push_back(triangle);		
				group->pushTriangle(triangle.fIndex);	

				while(iss >> str) {
					basis::StringLibrary::replace(str, "/", " ");
					iss_.clear(); iss_.str(str);
					if(!(iss_ >> vvi >> tti))
						break;

					triangle.vIndex[0] = triangles[triangle.fIndex].vIndex[0];
					triangle.tIndex[0] = triangles[triangle.fIndex].tIndex[0];
					triangle.vIndex[1] = triangles[triangle.fIndex].vIndex[2];
					triangle.tIndex[1] = triangles[triangle.fIndex].tIndex[2];
					triangle.vIndex[2] = vvi < 0 ? vvi + numvertices : vvi-1;
					triangle.tIndex[2] = tti < 0 ? tti + numtexcoords : tti-1;

					++triangle.fIndex;
					triangles.push_back(triangle);		
					group->pushTriangle(triangle.fIndex);	
				}
			}
			else {
				for(int i=0;i<3;++i) {	
					if(i != 0)
						iss >> str;

					iss_.clear(); iss_.str(str);
					iss_ >> vvi;
					triangle.vIndex[i] = vvi < 0 ? vvi + numvertices : vvi-1;
				}

				triangle.fIndex = static_cast<int>(triangles.size());
				triangles.push_back(triangle);		
				group->pushTriangle(triangle.fIndex);	

				while(iss >> str) {
					iss_.clear(); iss_.str(str);
					if(!(iss_ >> vvi))
						break;

					triangle.vIndex[0] = triangles[triangle.fIndex].vIndex[0];
					triangle.vIndex[1] = triangles[triangle.fIndex].vIndex[2];
					triangle.vIndex[2] = vvi < 0 ? vvi + numvertices : vvi-1;

					++triangle.fIndex;
					triangles.push_back(triangle);		
					group->pushTriangle(triangle.fIndex);	
				}
			}					
		}
	}

	numVertices = static_cast<int>(posVerts.size());
	numTriangles = static_cast<int>(triangles.size());
	numNormals = static_cast<int>(normals.size());
	numTexCoords = static_cast<int>(texCoords.size());
	numFaceNormals = static_cast<int>(faceNormals.size());
	numMaterials = static_cast<int>(materials.size());
	
}
