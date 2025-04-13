#ifndef _GEOMETRYSPACE_ADF_QUADTREE_H_
#define _GEOMETRYSPACE_ADF_QUADTREE_H_

#include "../BiTree/QuadTree_Vertices.h"
#include "../MarchingCube/MarchingTriangle.h"

namespace geometryspace {

template <typename T>
class ADFQuadTreeCell {
public : 
	ADFQuadTreeCell() : parent(NULL), children(NULL) {}	
	ADFQuadTreeCell(unsigned int _locCode[2], unsigned int _level, ADFQuadTreeCell* _parent, unsigned int vert[4])
	: level(_level), parent(_parent), children(NULL) { 
		memcpy(locCode, _locCode, 2 * sizeof(unsigned int));
		memcpy(vertices, vert, 4 * sizeof(unsigned int));
	}

public : 
	unsigned int		locCode[2];
	unsigned int		level;	
	ADFQuadTreeCell*	parent;	
	ADFQuadTreeCell*	children;

	unsigned int		vertices[4];
};

template <typename T>
class ADFQuadTree : public _QuadTree_Vertices<T, ADFQuadTreeCell<T>> {
public : 
	ADFQuadTree(unsigned int maxLevel = 2, Vector2<T> lc = Vector2<T>(0), Vector2<T> uc = Vector2<T>(1));

	void init();
	void init(unsigned int ml, Vector2<T> lc, Vector2<T> uc);

	template<typename DISTFIELD>
	void setPhi(const DISTFIELD &);

	//////////////////////////////////////////////////////////////////////////
	// Utility
	void marchingTriangle() const;

	//////////////////////////////////////////////////////////////////////////
	// glDraw
	void glDraw(const RGBA& color) const;

protected : 
	//////////////////////////////////////////////////////////////////////////
	// Basic Split Procedure
	template<typename _QUADTREE>
	class SplitProc : public _QuadTree_Vertices<T, ADFQuadTreeCell<T>>::SplitProc<_QUADTREE> {
	public : 
		SplitProc(_QUADTREE *qt) : _QuadTree_Vertices<T, ADFQuadTreeCell<T>>::SplitProc<_QUADTREE>(qt) {}

		ADFQuadTreeCell<T> * createRootCell(unsigned int rootLevel) const;

		inline const bool isSplit(ADFQuadTreeCell<T> *cell) const;

		void splitCell(ADFQuadTreeCell<T> *cell) const;
		void postProcess() const;
	};

	template<typename _QUADTREE, typename DISTFIELD> 
	class SplitProc_DistField : public SplitProc<_QUADTREE> {
	public : 
		SplitProc_DistField(_QUADTREE *qt, const DISTFIELD &_field) : SplitProc<_QUADTREE>(qt), field(_field) {}

		ADFQuadTreeCell<T> * createRootCell(unsigned int rootLevel) const;

		inline const bool isSplit(ADFQuadTreeCell<T> *cell) const;

		void splitCell(ADFQuadTreeCell<T> *cell) const;
		void postProcess() const;

	private : 
		const DISTFIELD &field;
	};

public : 
	typename ADFQuadTree<T>::SplitProc<ADFQuadTree> * getSplitProc() const;
	

public : 
	std::vector<T>	phi;

protected : 
	//////////////////////////////////////////////////////////////////////////
	// Variable for Visualization	
	mutable std::vector<Triangle2D<T>> renderingTriangles;
};

template <typename T>
ADFQuadTree<T>::ADFQuadTree(unsigned int maxLevel, Vector2<T> _lc, Vector2<T> _uc) 
: _QuadTree_Vertices<T,ADFQuadTreeCell<T>>(maxLevel, _lc, _uc) 
{
	init();
}

template <typename T>
void ADFQuadTree<T>::init() {
	_QuadTree_Vertices<T,ADFQuadTreeCell<T>>::init();
	phi.clear();
	renderingTriangles.clear();
}

template <typename T>
void ADFQuadTree<T>::init(unsigned int ml, Vector2<T> lc, Vector2<T> uc) {
	_QuadTree_Vertices<T,ADFQuadTreeCell<T>>::init(ml,lc,uc);
	phi.clear();
	renderingTriangles.clear();
}

template <typename T> template<typename DISTFIELD>
void ADFQuadTree<T>::setPhi(const DISTFIELD &field) {
	SplitProc_DistField<ADFQuadTree<T>,DISTFIELD> splitproc(this, field);

	init();
	makeTree(splitproc);
	marchingTriangle();
}

template <typename T> 
void ADFQuadTree<T>::marchingTriangle() const {
	renderingTriangles.clear();

	for(std::size_t i=0;i<numCells;++i) {

		if(!QuadTreeLibrary::isLeaf(cells[i]))
			continue;

		std::size_t iFlagIndex = 0, iEdgeFlags = 0;
		unsigned int nodes[4];

		nodes[0] = cells[i]->vertices[0];	nodes[1] = cells[i]->vertices[1];
		nodes[2] = cells[i]->vertices[3];	nodes[3] = cells[i]->vertices[2];

		for(std::size_t iVertex=0;iVertex<4;++iVertex) {				
			if(phi[nodes[iVertex]] <= 0) {
				iFlagIndex |= (1i64) << iVertex; 
			}
		}
		
		if(iFlagIndex == 15)	{
			Triangle2D<T> t;
			t.v[0] = posVerts[nodes[0]];	t.v[1] = posVerts[nodes[1]];	t.v[2] = posVerts[nodes[3]];	renderingTriangles.push_back(t);
			t.v[0] = posVerts[nodes[2]];	t.v[1] = posVerts[nodes[3]];	t.v[2] = posVerts[nodes[1]];	renderingTriangles.push_back(t);
		}
		if(iFlagIndex==0 || iFlagIndex==15) 
			continue;

		bool isVertexOnEdge[4];
		for(std::size_t l=0;l<4;l++)	isVertexOnEdge[l] = false;
		for(std::size_t l=0;EdgeConnectionTableMC[iFlagIndex][l]!=-1;l++)	isVertexOnEdge[EdgeConnectionTableMC[iFlagIndex][l]] = true;

		Vector2<T> cellsize = getCellSize(cells[i]->level);
		Vector2<T> vertices[4];
		for( std::size_t iEdge=0; iEdge<4; ++iEdge ) {		// Find the point of size_tersection of the surface per each edge
			if( isVertexOnEdge[iEdge] ) {					// If there is an intersection on this edge
				std::size_t e0 = EdgeConnection[iEdge][0];
				std::size_t e1 = EdgeConnection[iEdge][1];

				T alpha = abs( phi[nodes[e0]] ) / abs( phi[nodes[e1]] - phi[nodes[e0]] );

				vertices[iEdge][0] = posVerts[nodes[0]][0] + cellsize[0] * ( VertexOffset[e0][0] + alpha * EdgeDirection[iEdge][0] );
				vertices[iEdge][1] = posVerts[nodes[0]][1] + cellsize[1] * ( VertexOffset[e0][1] + alpha * EdgeDirection[iEdge][1] );
			}				
		}

		for(std::size_t m=0;m<10;m+=3)	{
			if(RenderingTriagleConnectionTableMC[iFlagIndex][m] == -1)
				break;

			Triangle2D<T> t;
			for(std::size_t l=0;l<3;++l)	{
				if(RenderingTriagleConnectionTableMC[iFlagIndex][m+l] >= 4 )	{
					t.v[l] = posVerts[nodes[RenderingTriagleConnectionTableMC[iFlagIndex][m+l]-4]];
				}			
				else	{
					std::size_t iVertex = RenderingTriagleConnectionTableMC[iFlagIndex][m+l];
					t.v[l] = vertices[iVertex];
				}			
			}
			renderingTriangles.push_back(t);
		}
	}
}

template <typename T> 
void ADFQuadTree<T>::glDraw(const RGBA& color) const {
	glColor4f(color[0], color[1], color[2], color[3]);

	glBegin(GL_TRIANGLES);
	for(std::size_t i=0;i<renderingTriangles.size();i++)	{
		glVertex3f(renderingTriangles[i].v[0][0], renderingTriangles[i].v[0][1], 0.0f);
		glVertex3f(renderingTriangles[i].v[1][0], renderingTriangles[i].v[1][1], 0.0f);
		glVertex3f(renderingTriangles[i].v[2][0], renderingTriangles[i].v[2][1], 0.0f);
	}
	glEnd();
}

//////////////////////////////////////////////////////////////////////////
// Basic Split Procedure
template <typename T> template<typename _QUADTREE>
ADFQuadTreeCell<T> * ADFQuadTree<T>::SplitProc<_QUADTREE>::createRootCell(unsigned int rootLevel) const {
	return _QuadTree_Vertices<T, ADFQuadTreeCell<T>>::SplitProc<_QUADTREE>::createRootCell(rootLevel);	
}

template <typename T> template<typename _QUADTREE>
inline const bool ADFQuadTree<T>::SplitProc<_QUADTREE>::isSplit(ADFQuadTreeCell<T> *cell) const {
	if(cell->level == 0)
		return false;
	return (rand() / (double)RAND_MAX) < 0.8f;
}

template <typename T> template<typename _QUADTREE>
void ADFQuadTree<T>::SplitProc<_QUADTREE>::splitCell(ADFQuadTreeCell<T> *cell) const {
	_QuadTree_Vertices<T, ADFQuadTreeCell<T>>::SplitProc<_QUADTREE>::splitCell(cell);
}

template <typename T> template<typename _QUADTREE>
void ADFQuadTree<T>::SplitProc<_QUADTREE>::postProcess() const {
	_QuadTree_Vertices<T, ADFQuadTreeCell<T>>::SplitProc<_QUADTREE>::postProcess();
}

template <typename T>
typename ADFQuadTree<T>::SplitProc<ADFQuadTree<T>> * ADFQuadTree<T>::getSplitProc() const {
	return new ADFQuadTree<T>::SplitProc<ADFQuadTree<T>>(const_cast<ADFQuadTree<T>*>(this));
}

//////////////////////////////////////////////////////////////////////////
// Split Procedure _ for Distance Fields
template <typename T> template<typename _QUADTREE, typename DISTFIELD>
ADFQuadTreeCell<T> * ADFQuadTree<T>::SplitProc_DistField<_QUADTREE,DISTFIELD>::createRootCell(unsigned int rootLevel) const {
	ADFQuadTreeCell<T> * root = _QuadTree_Vertices<T, ADFQuadTreeCell<T>>::SplitProc<_QUADTREE>::createRootCell(rootLevel);

	qt->phi.resize(4);
	for(size_t i=0;i<4;++i) 
		qt->phi[i] = field(qt->posVerts[i]);

	return root;
}

template <typename T> template<typename _QUADTREE, typename DISTFIELD>
inline const bool ADFQuadTree<T>::SplitProc_DistField<_QUADTREE,DISTFIELD>::isSplit(ADFQuadTreeCell<T> *cell) const {
	if(cell->level == 0)
		return false;

	Vector2<T> cellsize = qt->getCellSize(cell->level);
	T dist = field(qt->posVerts[cell->vertices[0]] + cellsize * T(0.5f));

	return (abs(dist) < cellsize[0]);
}

template <typename T> template<typename _QUADTREE, typename DISTFIELD>
void ADFQuadTree<T>::SplitProc_DistField<_QUADTREE,DISTFIELD>::splitCell(ADFQuadTreeCell<T> *cell) const {
	size_t preSize = qt->posVerts.size();
	_QuadTree_Vertices<T, ADFQuadTreeCell<T>>::SplitProc<_QUADTREE>::splitCell(cell);
	size_t afterSize = qt->posVerts.size();

	for(size_t i=preSize;i<afterSize;++i) 
		qt->phi.push_back(field(qt->posVerts[i]));
}

template <typename T> template<typename _QUADTREE, typename DISTFIELD>
void ADFQuadTree<T>::SplitProc_DistField<_QUADTREE,DISTFIELD>::postProcess() const {
	_QuadTree_Vertices<T, ADFQuadTreeCell<T>>::SplitProc<_QUADTREE>::postProcess();

//	cout << qt->posVerts.size() << endl;
//	cout << qt->phi.size() << endl;
}

};		// geometryspace

#endif