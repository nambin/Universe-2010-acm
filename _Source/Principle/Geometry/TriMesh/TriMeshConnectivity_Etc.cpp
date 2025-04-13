#include <stdafx.h>
#include <iostream>
#include <Principle/Geometry/Surface/TriMeshConnectivity.h>

using namespace geometryspace;

void TriMeshConnectivity::printV2F(int face) const {
	assert(face < numFaces);
	std::cout << "Face " << face << "th's Vertices : " 
			<< v2f(face,0) << " " << v2f(face,1) << " " << v2f(face,2) << std::endl;
}

bool TriMeshConnectivity::isClosedTwoManifold() const {
	for(int i=0;i<numEdges;++i) {
		if(f2e(i,0) == -1 || f2e(i,1) == -1)
			return false;
	}
	for(int i=0;i<numFaces;++i) {
		if(f2f(i,0) == -1 || f2f(i,1) == -1 || f2f(i,2) == -1)
			return false;
	}
	return true;
}