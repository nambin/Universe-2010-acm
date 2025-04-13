#include <stdafx.h>
#include <Principle/Geometry/Triangle/TriangleManipulator.h>

using namespace std;
using namespace geometryspace;

TriangleManipulator::TriangleManipulator() {
	has_Anglelimit = false;			 anglelimit = 0;
	has_Arealimit = false;			 arealimit = 0;
	encloseConvexHull = true;
	conformingDelaunay = false;
	edge_Output = true;
	voronoi_Output = false;
	neighbor_Output = true;
	prohibit_Steiner_all = false;		
	prohibit_Steiner_bdry = false;
	has_Steinerlimit = false;		steinerlimit = 0;
	quite_Message = true;
}

const string TriangleManipulator::getCommandLine() {
	char str[256];
	string command;

	command += "pz";
	
	if(has_Anglelimit) {
		if(anglelimit < 0)	command += "q";
		else {
			sprintf_s(str, "q%f", anglelimit);
			command += str;
		}
	}

	if(has_Arealimit) {
		if(arealimit < 0)	command += "a";
		else {
			sprintf_s(str, "a%f", arealimit);
			command += str;
		}
	}

	if(encloseConvexHull)	command += "c";
	if(conformingDelaunay)	command += "D";
	if(edge_Output)			command += "e";
	if(voronoi_Output)		command += "v";
	if(neighbor_Output)		command += "n";

	if(prohibit_Steiner_all)		command += "YY";
	else if(prohibit_Steiner_bdry)	command += "Y";

	if(has_Steinerlimit) {
		sprintf_s(str, "S%d", steinerlimit);
		command += str;
	}

	if(quite_Message)		command += "Q";

	return command;
}

void TriangleManipulator::makeTriangleMesh() {
	string command = getCommandLine();	

	char _command[1024];
	strcpy_s(_command, 1024, command.c_str());
	cout << "Command : " << _command << endl;

	// Define input
	struct triangulateio in;	

	// Points Initialization
	in.numberofpoints = static_cast<int>(vertices.size());
	in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
	for(int i=0;i<in.numberofpoints;++i) {
		in.pointlist[i*2+0] = vertices[i][0];
		in.pointlist[i*2+1] = vertices[i][1];
	}

	// Point Attributes Initialization
	in.numberofpointattributes = 0;
	in.pointattributelist = NULL;
	
	// Point Marker Initialization
	in.pointmarkerlist = (int *) malloc(in.numberofpoints * sizeof(int));
	for(int i=0;i<in.numberofpoints;++i)
		in.pointmarkerlist[i] = marker_vertices[i];

	in.numberofsegments = 0;
	in.numberofholes = 0;
	in.numberofregions = 0;
	in.regionlist = NULL;

	// Make necessary initializations so that 
	// Triangle can return a triangulation in `output' and 
	// a Voronoi Diagram in `vorout'.  
	struct triangulateio *output = new struct triangulateio();
	struct triangulateio *vorout = new struct triangulateio();

	output->pointlist = (REAL *) NULL;				// Not needed if -N switch used. 
	output->pointattributelist = (REAL *) NULL;		// Not needed if -N switch used or number of point attributes is zero: 
	output->pointmarkerlist = (int *) NULL;			// Not needed if -N or -B switch used. 
	output->trianglelist = (int *) NULL;			// Not needed if -E switch used. 
	output->triangleattributelist = (REAL *) NULL;	// Not needed if -E switch used or number of triangle attributes is zero: 
	output->neighborlist = (int *) NULL;			// Needed only if -n switch used. 
	output->segmentlist = (int *) NULL;				// Needed only if segments are output (-p or -c) and -P not used: 
	output->segmentmarkerlist = (int *) NULL;		// Needed only if segments are output (-p or -c) and -P and -B not used: 
	output->edgelist = (int *) NULL;				// Needed only if -e switch used. 
	output->edgemarkerlist = (int *) NULL;			// Needed if -e used and -B not used. 

	// Triangulate
	triangulate(_command, &in, output, vorout);

	// Copy
	TriangleManipulator::setMesh(mesh, *output);

	delete output;
	delete vorout;

	cout << "Number of Vertices : " << mesh.numVertices << endl;
	cout << "Number of Faces    : " << mesh.numFaces << endl;
}

void TriangleManipulator::makeRectangleMesh(const Vector2<double>& lb, const Vector2<double> ub) {
	vertices.push_back(lb);
	vertices.push_back(Vector2<double>(ub[0], lb[1]));	
	vertices.push_back(ub);
	vertices.push_back(Vector2<double>(lb[0], ub[1]));
	marker_vertices.push_back(1);
	marker_vertices.push_back(1);
	marker_vertices.push_back(1);
	marker_vertices.push_back(1);

	makeTriangleMesh();
}