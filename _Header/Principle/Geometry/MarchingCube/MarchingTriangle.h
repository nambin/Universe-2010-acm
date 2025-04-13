#ifndef _GEOMETRYSPACE_MARCHINGTRIANGLE_H_
#define _GEOMETRYSPACE_MARCHINGTRIANGLE_H_

#include "../../Arithmetic/Vector/Vector2.h"

namespace geometryspace {

using arithmetic::Vector2;

template<typename T>
class Triangle2D { 
public : 
	Vector2<T> v[3]; 
};

template<typename T>
class Segment2D {
public : 
	Segment2D() {}
	Segment2D(const Vector2<T>& v0, const Vector2<T>& v1) { v[0] = v0; v[1] = v1; }

	Vector2<T> v[2]; 
};

// VertexOffset lists the positions, relative to vertex0, of each of the 8 vertices of a cube
static const double VertexOffset[4][2] = {
	{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}
};

// EdgeConnection lists the index of the endpoint vertices for each of the 4 edges of the cube
static const int EdgeConnection[4][2] = { 
	{0,1}, {1,2}, {3,2}, {0,3}
};

// EdgeDirection lists the direction vector (vertex[i+1]-vertex[i])/2 for each edge in the cube
static const double EdgeDirection[4][2] = { 
	{1.0, 0.0}, {0.0, 1.0}, {1.0, 0.0}, {0.0, 1.0}
};

static const int EdgeConnectionTableMC[16][5] =	{
	{-1,-1,-1,-1,-1},	// 0
	{ 3, 0,-1,-1,-1},   // 1
	{ 0, 1,-1,-1,-1},   // 2
	{ 3, 1,-1,-1,-1},   // 3
	{ 1, 2,-1,-1,-1},	// 4
	{ 1, 2, 3, 0,-1},   // 5
	{ 0, 2,-1,-1,-1},   // 6
	{ 3, 2,-1,-1,-1},   // 7
	{ 2, 3,-1,-1,-1},   // 8
	{ 2, 0,-1,-1,-1},   // 9
	{ 2, 3, 0, 1,-1},   // 10
	{ 2, 1,-1,-1,-1},   // 11
	{ 1, 3,-1,-1,-1},   // 12
	{ 1, 0,-1,-1,-1},   // 13
	{ 0, 3,-1,-1,-1},   // 14
	{-1,-1,-1,-1,-1},	// 15
};

static const int RenderingTriagleConnectionTableMC[16][10] = {
	{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
	{ 3, 4, 0,-1,-1,-1,-1,-1,-1,-1},	// 1
	{ 0, 5, 1,-1,-1,-1,-1,-1,-1,-1},	// 2
	{ 3, 4, 1, 4, 5, 1,-1,-1,-1,-1},	// 3
	{ 1, 6, 2,-1,-1,-1,-1,-1,-1,-1},	// 4
	{ 1, 6, 2, 3, 4, 0,-1,-1,-1,-1},	// 5
	{ 0, 5, 2, 2, 5, 6,-1,-1,-1,-1},	// 6
	{ 3, 5, 2, 2, 5, 6, 3, 4, 5,-1},	// 7
	{ 2, 7, 3,-1,-1,-1,-1,-1,-1,-1},	// 8
	{ 2, 4, 0, 2, 7, 4,-1,-1,-1,-1},	// 9
	{ 2, 7, 3, 0, 5, 1,-1,-1,-1,-1},	// 10
	{ 2, 4, 1, 1, 4, 5, 2, 7, 4,-1},	// 11
	{ 1, 6, 3, 3, 6, 7,-1,-1,-1,-1},	// 12
	{ 1, 7, 0, 0, 7, 4, 1, 6, 7,-1},	// 13
	{ 0, 6, 3, 0, 5, 6, 3, 6, 7,-1},	// 14
	{ 4, 5, 6, 4, 6, 7,-1,-1,-1,-1},	// 15
};

template<typename T>
void contouring_bdry(std::vector<Segment2D<T>>& segments, T phi0, T phi1, const Vector2<T>& v0, const Vector2<T>& v1, const Vector2<T>& v_mid) {
	if(phi0 > 0 && phi1 > 0)	return;
	else if(phi0 > 0)			segments.push_back(Segment2D<T>(v_mid, v1));
	else if(phi1 > 0)			segments.push_back(Segment2D<T>(v_mid, v0));
	else						segments.push_back(Segment2D<T>(v0   , v1));
}

template<typename T>
void contouring(std::vector<Segment2D<T>>& segments, const T phi[4], const Vector2<T>& org, const Vector2<T>& h, bool bdry_xm, bool bdry_xp, bool bdry_ym, bool bdry_yp) {

	segments.clear();

	int iFlagIndex = 0, iEdgeFlags = 0;	

	for(int iVertex=0;iVertex<4;++iVertex) {				
		if(phi[iVertex] <= 0) {
			iFlagIndex |= (1i64) << iVertex; 
		}
	}

	if(iFlagIndex==0) 
		return;

	bool isVertexOnEdge[4];
	for(int l=0;l<4;l++)										isVertexOnEdge[l] = false;
	for(int l=0;EdgeConnectionTableMC[iFlagIndex][l]!=-1;l++)	isVertexOnEdge[EdgeConnectionTableMC[iFlagIndex][l]] = true;

	Vector2<T> vertices[4];
	for( int iEdge=0; iEdge<4; ++iEdge ) {				// Find the point of intersection of the surface per each edge
		if( isVertexOnEdge[iEdge] ) {					// If there is an intersection on this edge
			const int e0 = EdgeConnection[iEdge][0];
			const int e1 = EdgeConnection[iEdge][1];

			const T alpha = abs( phi[e0] ) / abs( phi[e1] - phi[e0] );

			vertices[iEdge][0] = org[0] + h[0] * ( VertexOffset[e0][0] + alpha * EdgeDirection[iEdge][0] );
			vertices[iEdge][1] = org[1] + h[1] * ( VertexOffset[e0][1] + alpha * EdgeDirection[iEdge][1] );
		}				
	}

	for(int m=0;m<5;m+=2)	{
		if(EdgeConnectionTableMC[iFlagIndex][m] == -1)
			break;

		Segment2D<T> seg;
		for(int l=0;l<2;++l) {
			const int iVertex = EdgeConnectionTableMC[iFlagIndex][m+l];
			seg.v[l] = vertices[iVertex];
		}
		segments.push_back(seg);
	}
	
	if(bdry_ym) { contouring_bdry(segments, phi[0], phi[1], org                        , org + Vector2<T>(h[0],  0 ), vertices[0]);  }
	if(bdry_xp) { contouring_bdry(segments, phi[1], phi[2], org + Vector2<T>(h[0],  0 ), org + Vector2<T>(h[0],h[1]), vertices[1]);  }
	if(bdry_yp) { contouring_bdry(segments, phi[2], phi[3], org + Vector2<T>(h[0],h[1]), org + Vector2<T>(  0 ,h[1]), vertices[2]);  }
	if(bdry_xm) { contouring_bdry(segments, phi[3], phi[0], org + Vector2<T>(  0 ,h[1]), org                        , vertices[3]);  }
}

template<typename T>
void contouring(std::vector<Triangle2D<T>>& triangles, const T phi[4], const Vector2<T>& org, const Vector2<T>& h) {

	triangles.clear();

	int iFlagIndex = 0, iEdgeFlags = 0;	

	for(int iVertex=0;iVertex<4;++iVertex) {				
		if(phi[iVertex] <= 0) {
			iFlagIndex |= (1i64) << iVertex; 
		}
	}

	if(iFlagIndex == 15)	{
		Triangle2D<T> t;
		t.v[0] = org;	  t.v[1] = org + Vector2<T>(h[0], 0);	t.v[2] = org + Vector2<T>(0, h[1]);	triangles.push_back(t);
		t.v[0] = org + h; t.v[1] = org + Vector2<T>(0, h[1]);	t.v[2] = org + Vector2<T>(h[0], 0);	triangles.push_back(t);
	}
	if(iFlagIndex==0 || iFlagIndex==15) 
		return;

	bool isVertexOnEdge[4];
	for(int l=0;l<4;l++)										isVertexOnEdge[l] = false;
	for(int l=0;EdgeConnectionTableMC[iFlagIndex][l]!=-1;l++)	isVertexOnEdge[EdgeConnectionTableMC[iFlagIndex][l]] = true;

	Vector2<T> vertices[4];
	for( int iEdge=0; iEdge<4; ++iEdge ) {		// Find the point of intersection of the surface per each edge
		if( isVertexOnEdge[iEdge] ) {					// If there is an intersection on this edge
			const int e0 = EdgeConnection[iEdge][0];
			const int e1 = EdgeConnection[iEdge][1];

			const T alpha = abs( phi[e0] ) / abs( phi[e1] - phi[e0] );

			vertices[iEdge][0] = org[0] + h[0] * ( VertexOffset[e0][0] + alpha * EdgeDirection[iEdge][0] );
			vertices[iEdge][1] = org[1] + h[1] * ( VertexOffset[e0][1] + alpha * EdgeDirection[iEdge][1] );
		}				
	}

	for(int m=0;m<10;m+=3)	{
		if(RenderingTriagleConnectionTableMC[iFlagIndex][m] == -1)
			break;

		Triangle2D<T> t;
		for(int l=0;l<3;++l)	{
			if(RenderingTriagleConnectionTableMC[iFlagIndex][m+l] >= 4 )	{
				t.v[l][0] = static_cast<T>(VertexOffset[RenderingTriagleConnectionTableMC[iFlagIndex][m+l]-4][0]) * h[0] + org[0];
				t.v[l][1] = static_cast<T>(VertexOffset[RenderingTriagleConnectionTableMC[iFlagIndex][m+l]-4][1]) * h[1] + org[1];
			}			
			else	{
				const int iVertex = RenderingTriagleConnectionTableMC[iFlagIndex][m+l];
				t.v[l] = vertices[iVertex];
			}			
		}
		triangles.push_back(t);
	}
}

};			// namespace geometryspace

#endif
