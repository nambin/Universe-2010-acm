#ifndef _GEOMETRYSPACE_TRIANGLEMESHVIEWER_H_
#define _GEOMETRYSPACE_TRIANGLEMESHVIEWER_H_

#include "../../Principle.h"
#include "../Triangle/TriangleManipulator.h"
#include "../../VirtualSpace/Scene/SceneObject.h"

namespace geometryspace {

using namespace basis;
using namespace arithmetic;
using namespace virtualspace;

class TriangleMeshViewer;
class PRINCIPLEDLL TriangleMeshViewer_IOHandler : public IOHandler {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	TriangleMeshViewer_IOHandler(TriangleMeshViewer *_viewer);
	~TriangleMeshViewer_IOHandler();

	//////////////////////////////////////////////////////////////////////////
	// Handler
	void keyboardHandler(unsigned char key, int	x, int y);
	void mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y);
	void mouseMotionHandler(int x, int y);

private : 
	TriangleMeshViewer *viewer;
};

class PRINCIPLEDLL TriangleMeshViewer : public SceneObject {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	TriangleMeshViewer();
	TriangleMeshViewer(const TriMesh<Vector2f> &_mesh);
	~TriangleMeshViewer();

	//////////////////////////////////////////////////////////////////////////
	// Setter
	void set(const TriMesh<Vector2f> &_mesh);

	//////////////////////////////////////////////////////////////////////////
	// glBind
	void glBind() const;

	//////////////////////////////////////////////////////////////////////////
	// Idle Function
	void idle();

	//////////////////////////////////////////////////////////////////////////
	// Initial View Setting
	const Vector3f viewEye() const { return Vector3f(0,0,2);	}
	const Vector3f viewRef() const { return Vector3f(0,0,0);	}
	const Vector3f viewUp() const  { return Vector3f(0,1,0);	}
	const float viewTheta() const  { return 45.0f;				}

public : 
	TriMesh<Vector2f> mesh;

	mutable unsigned int iv, nth;
};

};		// namespace geometryspace

#endif