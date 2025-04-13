#include <stdafx.h>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <vector>
#include <Principle/Geometry/Surface/TriMeshConnectivity.h>

using namespace std;
using namespace geometryspace;

// Set e2f, v2e, f2f
// Set f2e, f2v, e2v, v2v, isBdryVertex
void TriMeshConnectivity::makeConnectivityData() {

//	for(unsigned int i=0;i<v2f.size();++i) assert(v2f[i] < numVertices);

	makeEdgeData(e2f, v2e, f2f, v2f);

//	for(unsigned int i=0;i<e2f.size();++i) assert(e2f[i] < numEdges);
//	for(unsigned int i=0;i<v2e.size();++i) assert(v2e[i] < numVertices);

	makeF2E(f2e, e2f);
	makeVertexRelated(f2v, e2v, isBdryVertex, v2f, f2f);
	makeV2V(v2v, v2e);
}

// Set numEdges, e2f, v2e, f2f
void TriMeshConnectivity::makeEdgeData(Array2<int> &e2f, Array2<int> &v2e, Array2<int> &f2f, const Array2<int> &v2f) {

	Array<EdgeData> edgeList(numFaces*3);

	for(int i=0;i<numFaces;++i) 
		for(int k=0;k<3;++k) 
			edgeList[i*3+k].set(v2f(i,(k+1)%3), v2f(i,(k+2)%3), i, k);

	edgeList.sort();

	// Count Number of Edges
	int nEdges = 0;
	int v0 = -1, v1 = -1;
	for(int i=0;i<numFaces*3;++i) {
		if(v0 != edgeList[i].v0 || v1 != edgeList[i].v1) {
			v0 = edgeList[i].v0;
			v1 = edgeList[i].v1;
			++nEdges;
		}
	}

	// Initialization
	numEdges = nEdges;
	e2f.resize(numFaces,3);
	v2e.resize(numEdges,2);
	f2f.resize(numFaces,3);	f2f.set(-1);

	int idxEdge = -1, idxFace0, idxInFace0, idxFace1, idxInFace1;
	v0 = -1, v1 = -1;
	for(int i=0;i<numFaces*3;++i) {		
		if(v0 != edgeList[i].v0 || v1 != edgeList[i].v1) {
			v0 = edgeList[i].v0;
			v1 = edgeList[i].v1;
			++idxEdge;

			v2e[idxEdge*2+0] = v0;
			v2e[idxEdge*2+1] = v1;

			idxFace0 = edgeList[i].idxFace;
			idxInFace0 = edgeList[i].idxInFace;
		}
		else {
			idxFace1 = edgeList[i].idxFace;
			idxInFace1 = edgeList[i].idxInFace;

			f2f[idxFace0 * 3 + idxInFace0] = idxFace1;
			f2f[idxFace1 * 3 + idxInFace1] = idxFace0;
		}
		e2f[edgeList[i].idxFace * 3 + edgeList[i].idxInFace] = idxEdge;	
	}	
}

// Set f2e
void TriMeshConnectivity::makeF2E(Array2<int> &f2e, const Array2<int> &e2f) {

	Array<int> idx(numEdges); idx.set(0);
	f2e.resize(numEdges,2); f2e.set(-1);

	for(int i=0;i<numFaces;++i) {
		for(int k=0;k<3;++k) {
			f2e[e2f[i*3+k] * 2 + idx[e2f[i*3+k]]] = i;
			++idx[e2f[i*3+k]];

			// For Debug
			if(idx[e2f[i*3+k]] > 2) {
				int edge = e2f[i*3+k];
				cerr << f2e[edge*2] << " " << f2e[edge*2+1] << " " << i << endl;

				printV2F(f2e[edge*2]);
				printV2F(f2e[edge*2+1]);
				printV2F(i);
				system("pause");
			}
		}
	}

	// For Debug
	//for(int i=0;i<numEdges;++i)	{
	//	if(idx[i] != 2) {
	//		cerr << "Error in setf2e " << idx[i] << endl;
	//	}
	//}
}

// Set f2v, save cyclically
// Set isBdryVertex
// Set e2v, save cyclically 
void TriMeshConnectivity::makeVertexRelated(Array<vector<int>> &f2v, Array<vector<int>> &e2v, Array<bool> &isBdryVertex, const Array2<int> &v2f, const Array2<int> &f2f) {

	f2v.resize(numVertices);
	e2v.resize(numVertices);
	isBdryVertex.resize(numVertices);
	Array<vector<int>> f2vTmp(numVertices);

	for(int i=0;i<numFaces;++i) for(int k=0;k<3;++k) f2vTmp[v2f(i,k)].push_back(i);
	
	for(int i=0;i<numVertices;++i) {
		int iF_first = f2vTmp[i][0];
		int iF = iF_first;
		int iF_pre = -2;

		while(true) {			
			int iV_Local = -1;
			for(int k=0;k<3;++k) { if(v2f(iF,k) == i) { iV_Local = k; break; } }

			int iF01_Local = (iV_Local+1) % 3;	int iF01 = f2f(iF, iF01_Local);
			int iF02_Local = (iV_Local+2) % 3;	int iF02 = f2f(iF, iF02_Local);
			int iNeighborFace;

			if(iF_pre == -2) { if(iF01 != -1)	  { swap(iF01_Local, iF02_Local); swap(iF01, iF02); }}
			else			 { if(iF01 != iF_pre) { swap(iF01_Local, iF02_Local); swap(iF01, iF02); }}			

			iNeighborFace = iF02;

			if(iF01 != iF_pre && iF01 != -1)	iNeighborFace = iF01;
			else								iNeighborFace = iF02;

			if(iNeighborFace == -1)				{ isBdryVertex[i] = true;	break;	}
			else if(iNeighborFace == iF_first)	{ isBdryVertex[i] = false;	break;	}		
			else								{ iF_pre = iF; iF = iNeighborFace;	}
		}

		iF_first = iF;
		iF_pre = -2;			
		while(true) {
			int iV_Local = -1;
			for(int k=0;k<3;++k) { if(v2f(iF,k) == i) { iV_Local = k; break; } }

			int iF01_Local = (iV_Local+1) % 3;	int iF01 = f2f(iF, iF01_Local);
			int iF02_Local = (iV_Local+2) % 3;	int iF02 = f2f(iF, iF02_Local);
			int iNeighborFace_Local, iNeighborFace;

			if(iF_pre == -2) { if(iF01 != -1)	  { swap(iF01_Local, iF02_Local); swap(iF01, iF02); }}
			else			 { if(iF01 != iF_pre) { swap(iF01_Local, iF02_Local); swap(iF01, iF02); }}		
			
			iNeighborFace_Local = iF02_Local;
			iNeighborFace = iF02;

			if(iF_pre == -2)
				e2v[i].push_back(e2f(iF, iF01_Local));

			f2v[i].push_back(iF);

			if(iNeighborFace == -1)				{ e2v[i].push_back(e2f(iF,iNeighborFace_Local)); break;	}
			else if(iNeighborFace == iF_first) 	{ break;													}
			else								{ e2v[i].push_back(e2f(iF,iNeighborFace_Local)); iF_pre = iF; iF = iNeighborFace; }
		}
	}
}

void TriMeshConnectivity::makeV2V(Array<vector<int>> &v2v, const Array2<int> &v2e) {

	// Set v2v, Save cylinderically using v2e
	v2v.resize(numVertices);

	for(int i=0;i<numVertices;++i) {
		for(size_t j=0;j<e2v[i].size();++j) {
			int iV_neigh, iEdge = e2v[i][j];
			for(int k=0;k<2;++k) 
				if(i != v2e[iEdge*2+k])
					iV_neigh = v2e[iEdge*2+k];
			v2v[i].push_back(iV_neigh);
		}
	}

//	cout << "v2v[0] : " << v2v[0].size() << endl;
//	cout << "e2v[0] : " << e2v[0].size() << endl;
}
