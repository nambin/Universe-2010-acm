#ifndef _GEOMETRYSPACE_OCTREEVIEWER_H_
#define _GEOMETRYSPACE_OCTREEVIEWER_H_

#include "../../Principle.h"
#include "../../Arithmetic/Vector/Vector3.h"
#include "../../VirtualSpace/Scene/SceneObject.h"
#include "../DistanceFieldLibrary/ADFOctree.h"

namespace geometryspace {

using namespace basis;
using namespace arithmetic;
using namespace virtualspace;

class OctreeViewer;
class PRINCIPLEDLL OctreeViewer_IOHandler : public IOHandler {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	OctreeViewer_IOHandler(OctreeViewer *_viewer);
	~OctreeViewer_IOHandler();

	//////////////////////////////////////////////////////////////////////////
	// Handler
	void keyboardHandler(unsigned char key, int	x, int y);
	void mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y);
	void mouseMotionHandler(int x, int y);

private : 
	OctreeViewer *viewer;
};

class PRINCIPLEDLL OctreeViewer : public SceneObject {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	OctreeViewer();

	//////////////////////////////////////////////////////////////////////////
	// glBind
	void glBind() const;

	//////////////////////////////////////////////////////////////////////////
	// Idle Function
	void idle();

	//////////////////////////////////////////////////////////////////////////
	// Initial View Setting
	const Vector3f viewEye() const { return Vector3f(0.5f,0.5f,2.3f);	}
	const Vector3f viewRef() const { return Vector3f(0.5f,0.5f,0);		}
	const Vector3f viewUp() const  { return Vector3f(0,1,0);			}
	const float viewTheta() const  { return 45.0f;						}

public : 
	//Octree<double> octree;
	//Octree_Vertices<double> octree;
	ADFOctree<double> octree;

	Vector3<double> pos[100];
	unsigned int idxCell;

	bool drawTree;
	bool drawOnlyMinimumLeaf;
	bool drawRenderingPrimitive;
};

};		// namespace geometryspace

#endif