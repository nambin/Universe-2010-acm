#ifndef _GEOMETRYSPACE_TRIMESH_CONNECTIVITY_H_
#define _GEOMETRYSPACE_TRIMESH_CONNECTIVITY_H_

#include "../../../Basis/Basis/BasicLibrary.h"
#include "../../../Basis/Basis/Array.h"
#include "../../../Basis/Basis/Array2.h"
#include "../../Principle.h"

namespace geometryspace {

using std::vector;
using basis::Array;
using basis::Array2;
using basis::MIN;
using basis::MAX;

class PRINCIPLEDLL TriMeshConnectivity {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	TriMeshConnectivity();	
	TriMeshConnectivity(const TriMeshConnectivity *mesh);
	TriMeshConnectivity(const TriMeshConnectivity &mesh);
	virtual ~TriMeshConnectivity() {}

	//////////////////////////////////////////////////////////////////////////
	// Initialization
private :
	void init();
	void init(int nV, int nF, int nE);

public : 
	void set(const TriMeshConnectivity *mesh);
	void set(const TriMeshConnectivity &mesh);
	void set(int nV, const vector<int> &_v2f);
	void set(int nV, int nF, const int *_v2f);
	void set(int nV, const Array2<int> &_v2f);
	
	void release();	

	//////////////////////////////////////////////////////////////////////////
	//	Connectivity (Works for 2D-Manifolds)
	void	makeConnectivityData();
	
	void	makeEdgeData(Array2<int> &e2f, Array2<int> &v2e, Array2<int> &f2f, const Array2<int> &v2f);
	void	makeF2E(Array2<int> &f2e, const Array2<int> &e2f);
	void	makeVertexRelated(Array<vector<int>> &f2v, Array<vector<int>> &e2v, Array<bool> &isBdryVertex, const Array2<int> &v2f, const Array2<int> &f2f);	 // Works for 2D-Manifolds
	void	makeV2V(Array<vector<int>> &v2v, const Array2<int> &v2e);

	//////////////////////////////////////////////////////////////////////////
	// IO
	void	printV2F(int face) const;

	//////////////////////////////////////////////////////////////////////////
	// Manipulation
	bool	isClosedTwoManifold() const;

public : 
	int					numVertices;
	int					numEdges;
	int					numFaces;

	Array<bool>			isBdryVertex;	// [numVertices]

	Array2<int>			v2f;			// [numFaces*3]		vertices of each face	
	Array2<int>			e2f;			// [numFaces*3]		edges of each face
	Array2<int>			f2f;			// [numFaces*3]		neighbor faces to each face

	Array2<int>			v2e;			// [numEdges*2]		vertices of each edge
	Array2<int>			f2e;			// [numEdges*2]		neighbor faces to each edge				// Need to modify

	Array<vector<int>>	v2v;			// [?*numVertices]	neighbor vertices to each vertex
	Array<vector<int>>	e2v;			// [?*numVertices]	neighbor edges to each vertex
	Array<vector<int>>	f2v;			// [?*numVertices]	neighbor faces to each vertex	
};

struct EdgeData {
public : 
	EdgeData() {}
	EdgeData(int _v0, int _v1, int _idxFace, int _idxInFace) { set(_v0, _v1, _idxFace, _idxInFace);	}
	void set(int _v0, int _v1, int _idxFace, int _idxInFace) { v0 = MIN(_v0,_v1); v1 = MAX(_v0,_v1); idxFace = _idxFace; idxInFace = _idxInFace; }
	bool operator<(const EdgeData &e) const {
		if(v0 == e.v0)	return v1 < e.v1;
		else			return v0 < e.v0;
	}

public : 
	int v0, v1, idxFace, idxInFace;
};

};		// namespace geometryspace

#endif
