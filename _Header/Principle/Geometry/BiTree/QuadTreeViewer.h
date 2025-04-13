#ifndef _GEOMETRYSPACE_QUADTREEVIEWER_H_
#define _GEOMETRYSPACE_QUADTREEVIEWER_H_

#include "../../Principle.h"
#include "../../Arithmetic/Vector/Vector3.h"
#include "../../VirtualSpace/Scene/SceneObject.h"
#include "QuadTree_Vertices.h"
#include "QuadTree_Neighbors.h"
#include "../DistanceFieldLibrary/ADFQuadTree.h"

namespace geometryspace {

using namespace basis;
using namespace arithmetic;
using namespace virtualspace;

class QuadTreeViewer;
class PRINCIPLEDLL QuadTreeViewer_IOHandler : public IOHandler {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	QuadTreeViewer_IOHandler(QuadTreeViewer *_viewer);
	~QuadTreeViewer_IOHandler();

	//////////////////////////////////////////////////////////////////////////
	// Handler
	void keyboardHandler(unsigned char key, int	x, int y);
	void mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y);
	void mouseMotionHandler(int x, int y);

private : 
	QuadTreeViewer *viewer;
};

class PRINCIPLEDLL QuadTreeViewer : public SceneObject {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	QuadTreeViewer();

	//////////////////////////////////////////////////////////////////////////
	// glBind
	void glBind() const;

	//////////////////////////////////////////////////////////////////////////
	// Idle Function
	void idle();

	//////////////////////////////////////////////////////////////////////////
	// Initial View Setting
	const Vector3f viewEye() const { return Vector3f(0.5f,0.5f,1.43f);	}
	const Vector3f viewRef() const { return Vector3f(0.5f,0.5f,0);		}
	const Vector3f viewUp() const  { return Vector3f(0,1,0);			}
	const float viewTheta() const  { return 45.0f;						}

public : 
	//QuadTree<double> quadtree;
	//QuadTree_Vertices<double> quadtree;
	//QuadTree_Neighbors<double> quadtree;
	ADFQuadTree<double> quadtree;

	Vector2<double> pos[100];
	unsigned int idxCell;

	bool drawTree;
	bool drawRenderingPrimitive;
};

};		// namespace geometryspace

#endif