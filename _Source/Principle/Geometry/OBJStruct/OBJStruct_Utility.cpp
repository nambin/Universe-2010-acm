#include <stdafx.h>
#include <Principle/Geometry/Surface/OBJStruct.h>

using namespace geometryspace;

OBJGroup* OBJStruct::addGroup(const char* const name) {

	OBJGroup* group = findGroup(name);
	if (group == NULL) {
		group = new OBJGroup(name);

		group->next = groups;
		groups = group;
		++numGroups;
	}

	return group;
}

OBJGroup* OBJStruct::findGroup(const char* const name) {

	OBJGroup* group = groups;
	while(group != NULL) {
		if(group->name.compare(name)==0)
			break;
		group = group->next;
	}

	return group;
}

const int OBJStruct::findMaterial(const char* const name)	{

	for (int i = 0; i < numMaterials; ++i) {
		if(materials[i].name.compare(name)==0)
			return i;
	}

	std::cout << "OBJStruct::findMaterial() : can't find material _ " << name << endl;
	return 0;		
}

void OBJStruct::setFaceNormals() { 

	faceNormals.clear();

	Vector3<double> n, v0, v1, v2;
	for(int i=0;i<numTriangles;++i) {
		const OBJTriangle &triangle = triangles[i];
		v0 = posVerts[triangle.vIndex[0]];
		v1 = posVerts[triangle.vIndex[1]];
		v2 = posVerts[triangle.vIndex[2]];

		n = arithmetic::unitVector(arithmetic::cross(v1-v0, v2-v0));
		faceNormals.push_back(Vector3<float>(n));
	}
}

void OBJStruct::unitize() {

	float maxx, minx, maxy, miny, maxz, minz;

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
	float w = abs(maxx-minx);
	float h = abs(maxy-miny);
	float d = abs(maxz-minz);

	// calculate center of the model
	float cx = (maxx + minx) * 0.5f;
	float cy = (maxy + miny) * 0.5f;
	float cz = (maxz + minz) * 0.5f;

	// calculate unitizing scale factor 
	float s = 1.0f / basis::MAX3(w, h, d);

	// translate around center then scale 
	translate(-Vector3<float>(cx, cy, cz));
	scale(s);	

	// translate around minimum vertex then scale 
//	translate(-Vector3<float>(minx, miny, minz));
//	scale(s);	
}

void OBJStruct::translate(Vector3<float> trans) {
	for(int i = 0; i < numVertices; ++i) 
		posVerts[i] += trans;
}

void OBJStruct::scale(float s)	{
	for(int i = 0; i < numVertices; ++i) {
		posVerts[i][0] *= s;
		posVerts[i][1] *= s;
		posVerts[i][2] *= s;
	}
}
