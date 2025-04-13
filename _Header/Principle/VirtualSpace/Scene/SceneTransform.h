#ifndef _VIRTUALSPACE_SCENETRANSFORM_H_
#define _VIRTUALSPACE_SCENETRANSFORM_H_

#include "SceneStateNode.h"
#include "../../Arithmetic/Transform/Transform.h"

namespace virtualspace {

using arithmetic::Transform;
using arithmetic::Matrix4f;

class Scene;
class PRINCIPLEDLL SceneTransform : public SceneStateNode {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	SceneTransform() : SceneStateNode(NULL)								{ }
	SceneTransform(const SceneTransform &src) : SceneStateNode(NULL)	{ setTransform(src);			}
	SceneTransform(const Matrix4f &src)	: SceneStateNode(NULL)			{ setTransform(src);			}
	SceneTransform(const Transform<float> &src) : SceneStateNode(NULL)	{ setTransform(src);			}
	~SceneTransform() {};

	//////////////////////////////////////////////////////////////////////////
	// Setter
	void setTransform(const SceneTransform &node)	{ transform = node.transform;	}	
	void setTransform(const Matrix4f &mat)			{ transform.set(mat);			}
	void setTransform(const Transform<float> &trans){ transform.set(trans);		}

	//////////////////////////////////////////////////////////////////////////
	// Getter
	const Transform<float> &getTransform() const	{ return transform;	}
	inline void glTransform() const					{ glMultMatrixf(transform.getMatrixPtr()->getptr()); }

	//////////////////////////////////////////////////////////////////////////
	// gl Operation
	void glBind() const;
	void glUnbind() const;

	//////////////////////////////////////////////////////////////////////////
	// Idle Function
	void idle();	

private : 
	Transform<float> transform;
};

};	// namespace virtualspace

#endif