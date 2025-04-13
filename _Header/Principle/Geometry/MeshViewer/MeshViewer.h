#ifndef _GEOMETRYSPACE_MESHVIEWER_H_
#define _GEOMETRYSPACE_MESHVIEWER_H_

#include "../../Principle.h"
#include "../Surface/TriMesh.h"
#include "../Surface/TriMesh_Simple.h"
#include "../Surface/OBJStruct.h"
#include "../../VirtualSpace/Scene/SceneObject.h"

namespace geometryspace {

using namespace basis;
using namespace arithmetic;
using namespace virtualspace;

class MeshViewer;
class PRINCIPLEDLL MeshViewer_IOHandler : public IOHandler {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	MeshViewer_IOHandler(MeshViewer *_viewer);
	~MeshViewer_IOHandler();

	//////////////////////////////////////////////////////////////////////////
	// Handler
	void keyboardHandler(unsigned char key, int	x, int y);
	void mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y);
	void mouseMotionHandler(int x, int y);

private : 
	MeshViewer *viewer;
};

class PRINCIPLEDLL MeshViewer : public SceneObject {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	MeshViewer();
	MeshViewer(const TriMesh<Vector3f> &_mesh);
	~MeshViewer();

	//////////////////////////////////////////////////////////////////////////
	// Setter
	void set(const TriMesh<Vector3f> &_mesh);

	//////////////////////////////////////////////////////////////////////////
	//	Loader
	void load(const char * const filename, bool change_order = false);

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
	TriMesh_Simple<Vector3f, vector<Vector3f>> mesh_simple;
	TriMesh<Vector3f> mesh;
	OBJStruct objStruct;
	PLYStruct plyStruct;

	int selected_face;
};

};		// namespace geometryspace

#endif