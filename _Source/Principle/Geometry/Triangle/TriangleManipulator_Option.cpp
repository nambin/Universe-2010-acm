#include <stdafx.h>
#include <Principle/Geometry/Triangle/TriangleManipulator.h>

using namespace std;
using namespace geometryspace;

void TriangleManipulator::setAnglelimit(double limit) {
	if(limit <= 0)	has_Anglelimit = false;
	else {
		has_Anglelimit = true;
		anglelimit = limit;
	}
}

void TriangleManipulator::setArealimit(double limit) {
	if(limit < 0)	has_Arealimit = false;
	else {
		has_Arealimit = true;
		arealimit = limit;
	}
}

void TriangleManipulator::setSteinerlimit(size_t limit) {
	if(limit == 0)	has_Steinerlimit = false;
	else {
		has_Steinerlimit = true;
		steinerlimit = limit;
	}
}

void TriangleManipulator::setFlag_ConvexHull(bool flag)					{ encloseConvexHull = flag;		}
void TriangleManipulator::setFlag_CDT(bool flag)						{ conformingDelaunay = flag;	}
void TriangleManipulator::setFlag_EdgeOutput(bool flag)					{ edge_Output = flag;			}
void TriangleManipulator::setFlag_VoronoiOutput(bool flag)				{ voronoi_Output = flag;		}
void TriangleManipulator::setFlag_NeighborOutput(bool flag)				{ neighbor_Output = flag;		}
void TriangleManipulator::setFlag_Prohibit_Steiner_Boundary(bool flag)	{ prohibit_Steiner_bdry = flag;	}
void TriangleManipulator::setFlag_Prohibit_Steiner_All(bool flag)		{ prohibit_Steiner_all = flag;	}
void TriangleManipulator::setFlag_Quite(bool flag)						{ quite_Message = flag;			}


