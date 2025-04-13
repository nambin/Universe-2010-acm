#ifndef _GEOMETRY_SPACE_TRIANGLE_MANIPULATOR_H_
#define _GEOMETRY_SPACE_TRIANGLE_MANIPULATOR_H_

#include <string>
#include <vector>
#include "../../Arithmetic/Vector/Vector2.h"
#include "triangle.h"
#include "../Surface/TriMesh.h"
#include "../../Principle.h"

namespace geometryspace {

using std::string;
using std::vector;
using arithmetic::Vector2;

extern "C" void triangulate(char *, struct triangulateio *, struct triangulateio *, struct triangulateio *);

// Alliez Optimization
// Point Attribute
// Triangle Attribute
// Hole
// Region
// Segments (for non convex)

// Refine

class PRINCIPLEDLL TriangleManipulator {
public : 
	TriangleManipulator();

	void setAnglelimit(double limit);
	void setArealimit(double limit);
	void setFlag_ConvexHull(bool);
	void setFlag_CDT(bool);
	void setFlag_EdgeOutput(bool);
	void setFlag_VoronoiOutput(bool);
	void setFlag_NeighborOutput(bool);
	void setFlag_Prohibit_Steiner_Boundary(bool);
	void setFlag_Prohibit_Steiner_All(bool);
	void setSteinerlimit(size_t limit);
	void setFlag_Quite(bool);

	const string getCommandLine();
	void makeTriangleMesh();

	void makeRectangleMesh(const Vector2<double>& lb, const Vector2<double> ub);

	template<typename T> static void setMesh(TriMesh<Vector2<T>>&, const triangulateio&);

public : 
	TriMesh<Vector2<double>>	mesh;

	vector<Vector2<double>>		vertices;
	vector<int>					marker_vertices;

	bool has_Anglelimit;		double anglelimit;
	bool has_Arealimit;			double arealimit;
	bool encloseConvexHull;
	bool conformingDelaunay;
	bool edge_Output;
	bool voronoi_Output;
	bool neighbor_Output;
	bool prohibit_Steiner_bdry;			
	bool prohibit_Steiner_all;
	bool has_Steinerlimit;		size_t steinerlimit;
	bool quite_Message;
};

#pragma warning(disable:4800)   

template<typename T> 
void TriangleManipulator::setMesh(TriMesh<Vector2<T>>& mesh, const triangulateio& tri) {
	mesh.release();

	// Assign Information of Vertices, Faces
	mesh.numVertices = tri.numberofpoints;
	mesh.numFaces = tri.numberoftriangles;

	mesh.posVerts.resize(mesh.numVertices);
	mesh.v2f.resize(mesh.numFaces, 3);

	for(int i=0;i<mesh.numVertices;++i)
		mesh.posVerts[i] = Vector2<T>(tri.pointlist[i*2], tri.pointlist[i*2+1]);

	for(int i=0;i<mesh.numFaces;++i) {
		mesh.v2f(i,0) = tri.trianglelist[i*3+0];
		mesh.v2f(i,1) = tri.trianglelist[i*3+1];
		mesh.v2f(i,2) = tri.trianglelist[i*3+2];
	}

	// Assign Marker Information of Vertices, if exists
	if(tri.pointmarkerlist != NULL) {
		mesh.isBdryVertex.resize(mesh.numVertices);
		for(int i=0;i<mesh.numVertices;++i)
			mesh.isBdryVertex[i] = static_cast<bool>(tri.pointmarkerlist[i]);
	}

	// Assign Information of Edge
	if(tri.edgelist != NULL) {
		mesh.numEdges = tri.numberofedges;
		mesh.v2e.resize(mesh.numEdges, 2);
		for(int i=0;i<mesh.numEdges;++i) {
			mesh.v2e(i,0) = tri.edgelist[i*2+0];
			mesh.v2e(i,1) = tri.edgelist[i*2+1];
		}
	}

	// Assign Neighborhood Information 
	if(tri.neighborlist != NULL) {
		mesh.f2f.resize(mesh.numFaces, 3);
		for(int i=0;i<mesh.numFaces;++i) {
			mesh.f2f(i,0) = tri.neighborlist[i*3+0];
			mesh.f2f(i,1) = tri.neighborlist[i*3+1];
			mesh.f2f(i,2) = tri.neighborlist[i*3+2];
		}
	}
}

};		// geometryspace

#endif