#ifndef _GEOMETRYSPACE_ADF_OCTREE_H_
#define _GEOMETRYSPACE_ADF_OCTREE_H_

#include <vector>
#include "../BiTree/Octree_Vertices.h"
#include "../MarchingCube/MarchingCube.h"
#include "../Surface/TriMesh_Simple.h"

namespace geometryspace {

using std::vector;

template <typename T>
class ADFOctreeCell {
public : 
	ADFOctreeCell() : parent(NULL), children(NULL) {}	
	ADFOctreeCell(unsigned int _locCode[3], unsigned int _level, ADFOctreeCell* _parent, unsigned int vert[8])
	: level(_level), parent(_parent), children(NULL) { 
		memcpy(locCode, _locCode, 3 * sizeof(unsigned int));
		memcpy(vertices, vert, 8 * sizeof(unsigned int));
	}

public : 
	unsigned int		locCode[3];
	unsigned int		level;	
	ADFOctreeCell*		parent;	
	ADFOctreeCell*		children;

	unsigned int		vertices[8];
};

template <typename T>
class ADFOctree : public _Octree_Vertices<T, ADFOctreeCell<T>> {
public : 
	ADFOctree(unsigned int maxLevel = 2, Vector3<T> _lc = Vector3<T>(0), Vector3<T> _uc = Vector3<T>(1));

	void init();
	void init(unsigned int ml, Vector3<T> lc, Vector3<T> uc);

	template<typename DISTFIELD>
	void setPhi(const DISTFIELD &);

	//////////////////////////////////////////////////////////////////////////
	// Utility
	void marchingCube() const;

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	void glDraw(const RGBA& color) const;

protected : 
	//////////////////////////////////////////////////////////////////////////
	// Basic Split Procedure
	template<typename _OCTREE>
	class SplitProc : public _Octree_Vertices<T, ADFOctreeCell<T>>::SplitProc<_OCTREE> {
	public : 
		SplitProc(_OCTREE *oct) : _Octree_Vertices<T, ADFOctreeCell<T>>::SplitProc<_OCTREE>(oct) {}

		ADFOctreeCell<T> * createRootCell(unsigned int rootLevel) const;

		inline const bool isSplit(ADFOctreeCell<T> *cell) const;

		void splitCell(ADFOctreeCell<T> *cell) const;
		void postProcess() const;
	};

	template<typename _OCTREE, typename DISTFIELD> 
	class SplitProc_DistField : public SplitProc<_OCTREE> {
	public : 
		SplitProc_DistField(_OCTREE *oct, const DISTFIELD &_field) : SplitProc<_OCTREE>(oct), field(_field) {}

		ADFOctreeCell<T> * createRootCell(unsigned int rootLevel) const;

		inline const bool isSplit(ADFOctreeCell<T> *cell) const;

		void splitCell(ADFOctreeCell<T> *cell) const;
		void postProcess() const;

	private : 
		const DISTFIELD &field;
	};

public : 
	typename ADFOctree<T>::SplitProc<ADFOctree> * getSplitProc() const;

public : 
	std::vector<T>	phi;

protected :
	//////////////////////////////////////////////////////////////////////////
	// Variable for Visualization	
	mutable TriMesh_Simple<Vector3<T>, vector<Vector3<T>>>		mesh;	
	mutable stdext::hash_map<unsigned long long, unsigned int>	edgeList;
};

template <typename T>
ADFOctree<T>::ADFOctree(unsigned int maxLevel, Vector3<T> lc, Vector3<T> uc) 
: _Octree_Vertices<T,ADFOctreeCell<T>>(maxLevel, lc, uc) 
{
	init();
}

template <typename T>
void ADFOctree<T>::init() {
	_Octree_Vertices<T,ADFOctreeCell<T>>::init();
	phi.clear();
	mesh.init();
	edgeList.clear();
}

template <typename T>
void ADFOctree<T>::init(unsigned int ml, Vector3<T> lc, Vector3<T> uc) {
	_Octree_Vertices<T,ADFOctreeCell<T>>::init(ml,lc,uc);
	phi.clear();
	mesh.init();
	edgeList.clear();
}

template <typename T> template<typename DISTFIELD>
void ADFOctree<T>::setPhi(const DISTFIELD &field) {
	SplitProc_DistField<ADFOctree<T>,DISTFIELD> splitproc(this, field);

	init();
	makeTree(splitproc);
	marchingCube();
}

template <typename T> 
void ADFOctree<T>::marchingCube() const {
	mesh.init();
	edgeList.clear();

	unsigned int nodes[8];
	Vector3<T> pt[12];
	unsigned long long ei;
	unsigned int newPtOnEdge[12];

	for(size_t i=0;i<cells.size();++i) {

		if(!OctreeLibrary::isLeaf(cells[i]))
			continue;

		int iFlagIndex = 0, iEdgeFlags = 0;

		nodes[0] = cells[i]->vertices[0]; nodes[1] = cells[i]->vertices[1]; nodes[2] = cells[i]->vertices[5]; nodes[3] = cells[i]->vertices[4];
		nodes[4] = cells[i]->vertices[2]; nodes[5] = cells[i]->vertices[3]; nodes[6] = cells[i]->vertices[7]; nodes[7] = cells[i]->vertices[6];

		for(int iVertex=0;iVertex<8;iVertex++) { if(phi[nodes[iVertex]]<=0) { iFlagIndex |= 1<<iVertex; } } // Which vertices are inside? If i-th vertex is inside, mark '1' at i-th bit of 'iFlagIndex'.

		if(iFlagIndex==0 || iFlagIndex==255) 
			continue; 

		iEdgeFlags = CubeEdgeFlags3D[iFlagIndex]; // Which edges intersect the surface? If i-th edge intersects the surface, mark '1' at i-th bit of 'iEdgeFlags' 
		
		Vector3<T> cellsize = getCellSize(cells[i]->level);

		T EPS = 1.e-3;
		for( int iEdge=0; iEdge<12; iEdge++ ) { // Find the point of intersection of the surface per each edge

			if( iEdgeFlags & (1<<iEdge) ) { // If there is an intersection on this edge

				int e0 = EdgeConnection3D[iEdge][0];
				int e1 = EdgeConnection3D[iEdge][1];

				T alpha = abs( phi[nodes[e0]] ) / abs( phi[nodes[e1]] - phi[nodes[e0]] );

				if(alpha < EPS)				alpha = EPS;
				else if(alpha > 1.0f - EPS)	alpha = 1.0f - EPS;

				pt[iEdge][0] = posVerts[nodes[0]][0] + cellsize[0] * ( VertexOffset3D[e0][0] + alpha * EdgeDirection3D[iEdge][0] );
				pt[iEdge][1] = posVerts[nodes[0]][1] + cellsize[1] * ( VertexOffset3D[e0][1] + alpha * EdgeDirection3D[iEdge][1] );
				pt[iEdge][2] = posVerts[nodes[0]][2] + cellsize[2] * ( VertexOffset3D[e0][2] + alpha * EdgeDirection3D[iEdge][2] );

				unsigned int a = nodes[e0], b = nodes[e1];
				if(a > b) swap(a, b);

				ei = a;
				ei = ei << 32;
				ei += b;

				if(edgeList.count(ei))
					newPtOnEdge[iEdge] = edgeList[ei];
				else {
					mesh.posVerts.push_back( pt[iEdge] );
					newPtOnEdge[iEdge] = mesh.numVertices++;
					edgeList[ei] = newPtOnEdge[iEdge];
				}
			}
		}

		int threeVerts[3];

		for( int iTriangle=0; iTriangle<5; ++iTriangle ) {

			// If there isn't any triangle to be drawn, escape this loop.
			if( TriangleConnectionTableMC3D[iFlagIndex][3*iTriangle] < 0 ) break;

			for( int iCorner=0; iCorner<3; iCorner++ ) {
				int iVertex = TriangleConnectionTableMC3D[iFlagIndex][3*iTriangle+iCorner];						
				threeVerts[iCorner] = newPtOnEdge[iVertex];
			}

			if(threeVerts[0] == -1 || threeVerts[1] == -1 || threeVerts[2] == -1)	{
				cout << "Error in MC reconstruction !" << endl;
				cout << i << " " << iFlagIndex << endl;
				system("pause");
			}

			mesh.v2f.push_back( threeVerts[0] );
			mesh.v2f.push_back( threeVerts[1] );
			mesh.v2f.push_back( threeVerts[2] );
			mesh.numFaces++;
		}
	}

	mesh.makeFaceNormals();
	mesh.makeVertexNormals();
}

template <typename T> 
void ADFOctree<T>::glDraw(const RGBA& color) const {
	mesh.glDraw_SmoothSurface(color);
}

//////////////////////////////////////////////////////////////////////////
// Basic Split Procedure
template <typename T> template<typename _OCTREE>
ADFOctreeCell<T> * ADFOctree<T>::SplitProc<_OCTREE>::createRootCell(unsigned int rootLevel) const {
	return _Octree_Vertices<T, ADFOctreeCell<T>>::SplitProc<_OCTREE>::createRootCell(rootLevel);	
}

template <typename T> template<typename _OCTREE>
inline const bool ADFOctree<T>::SplitProc<_OCTREE>::isSplit(ADFOctreeCell<T> *cell) const {
	if(cell->level == 0)
		return false;
	return (rand() / (double)RAND_MAX) < 0.4f;
}

template <typename T> template<typename _OCTREE>
void ADFOctree<T>::SplitProc<_OCTREE>::splitCell(ADFOctreeCell<T> *cell) const {
	_Octree_Vertices<T, ADFOctreeCell<T>>::SplitProc<_OCTREE>::splitCell(cell);
}

template <typename T> template<typename _OCTREE>
void ADFOctree<T>::SplitProc<_OCTREE>::postProcess() const {
	_Octree_Vertices<T, ADFOctreeCell<T>>::SplitProc<_OCTREE>::postProcess();
}

template <typename T>
typename ADFOctree<T>::SplitProc<ADFOctree<T>> * ADFOctree<T>::getSplitProc() const {
	return new ADFOctree<T>::SplitProc<ADFOctree<T>>(const_cast<ADFOctree<T>*>(this));
}

//////////////////////////////////////////////////////////////////////////
// Split Procedure _ for Distance Fields
template <typename T> template<typename _OCTREE, typename DISTFIELD>
ADFOctreeCell<T> * ADFOctree<T>::SplitProc_DistField<_OCTREE,DISTFIELD>::createRootCell(unsigned int rootLevel) const {
	ADFOctreeCell<T> * root = _Octree_Vertices<T, ADFOctreeCell<T>>::SplitProc<_OCTREE>::createRootCell(rootLevel);

	oct->phi.resize(8);
	for(size_t i=0;i<8;++i) 
		oct->phi[i] = field(oct->posVerts[i]);

	return root;
}

template <typename T> template<typename _OCTREE, typename DISTFIELD>
inline const bool ADFOctree<T>::SplitProc_DistField<_OCTREE,DISTFIELD>::isSplit(ADFOctreeCell<T> *cell) const {
	if(cell->level == 0)
		return false;

	Vector3<T> cellsize = oct->getCellSize(cell->level);
	T dist = field(oct->posVerts[cell->vertices[0]] + cellsize * T(0.5f));

	return (abs(dist) < cellsize[0]);
}

template <typename T> template<typename _OCTREE, typename DISTFIELD>
void ADFOctree<T>::SplitProc_DistField<_OCTREE,DISTFIELD>::splitCell(ADFOctreeCell<T> *cell) const {
	size_t preSize = oct->posVerts.size();
	_Octree_Vertices<T, ADFOctreeCell<T>>::SplitProc<_OCTREE>::splitCell(cell);
	size_t afterSize = oct->posVerts.size();

	for(size_t i=preSize;i<afterSize;++i) 
		oct->phi.push_back(field(oct->posVerts[i]));
}

template <typename T> template<typename _OCTREE, typename DISTFIELD>
void ADFOctree<T>::SplitProc_DistField<_OCTREE,DISTFIELD>::postProcess() const {
	_Octree_Vertices<T, ADFOctreeCell<T>>::SplitProc<_OCTREE>::postProcess();

	//	cout << oct->posVerts.size() << endl;
	//	cout << oct->phi.size() << endl;
}

};		// geometryspace

#endif