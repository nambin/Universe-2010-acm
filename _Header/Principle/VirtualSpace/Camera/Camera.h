#ifndef _VIRTUALSPACE_CAMERA_H_
#define _VIRTUALSPACE_CAMERA_H_

#include "../../Arithmetic/Vector/Vector3.h"
#include "../../Arithmetic/Orientation/Quaternion.h"
#include "../RenderingPipeline/RenderingPipelineSetting.h"
#include "../Scene/SceneObject.h"
#include "../../Principle.h"

namespace virtualspace {

using namespace basis;
using arithmetic::Vector3f;
using arithmetic::Matrix4f;

class Camera;
class PRINCIPLEDLL Camera_IOHandler : public IOHandler {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Camera_IOHandler(Camera *c);
	~Camera_IOHandler();

	//////////////////////////////////////////////////////////////////////////
	// Handler
	void keyboardHandler(unsigned char key, int	x, int y);
	void mouseHandler(IO_MOUSE_BUTTON button, IO_MOUSE_STATE state, int x, int y);
	void mouseMotionHandler(int x, int y);

private : 
	Camera *camera;
};

class PRINCIPLEDLL Camera : public SceneObject {
public : 
	enum CC_MODE { CC_IDLE_MODE = 0, CC_TRACKBALL_MODE, CC_TRANS_MODE, CC_PAN_MODE };

public : 
	struct Ray {
		Ray() {}
		Ray(const Vector3f &_org, const Vector3f &_dir) {
			org = _org; dir = _dir;
		}

		Vector3f org, dir;
	};

public:
	//////////////////////////////////////////////////////////////////////////
	// Constructor
	Camera();
	Camera(Vector3f _eye, Vector3f _ref, Vector3f _up, float _theta);
	void set();
	void set(Vector3f _eye, Vector3f _ref, Vector3f _up, float _theta);	

	virtual void glBind() const;
	void glDraw() const {}
	virtual void idle() {}

	//////////////////////////////////////////////////////////////////////////
	// GL rendering
	void			lookAt() const;
	void			viewport() const;
	void			perspective() const;
	void			resizeScreen(int _width, int _height);

	const Matrix4f	getModelViewMatrix() const;
	const int		getScreenWidth()	 const { return screenWidth;	}
	const int		getScreenHeight()	 const { return screenHeight;	}

	//////////////////////////////////////////////////////////////////////////
	// Trackball related
	void			trackballMouseDn(int x, int y);
	void			trackballMouseUp(int x, int y);
	void			trackballMouseMove(int x, int y);

	//////////////////////////////////////////////////////////////////////////
	// Translate related 
	void			translateMouseDn(int x, int y);
	void			translateMouseUp(int x, int y);
	void			translateMouseMove(int x, int y);

	//////////////////////////////////////////////////////////////////////////
	// Pan in/out
	void			pan(float mag);
	void			panin();
	void			panout();

	void			panMouseDn(int x, int y);
	void			panMouseUp(int x, int y);
	void			panMouseMove(int x, int y);

	//////////////////////////////////////////////////////////////////////////
	// General function
	const Vector3f	hemisphereMap(int x, int y) const;
	void			mouseMove(int x, int y);
	const Ray		getLastRay(int x, int y) const;
	const Ray		getRay(int x, int y) const ;

	static const float PanSpeed()		{ return 0.05f;		}
	static const float RotationSpeed()	{ return 1.0f;		}
	static const float TransSpeed()		{ return 0.0005f;	}
	static const float PerspectiveNear() { return 0.001f;	}
	static const float PerspectiveFar()	{ return 50000.0f;	}

	void printCamera() const;

	const Vector3<float> viewEye() const { return eye; }
	const Vector3<float> viewRef() const { return ref; }
	const Vector3<float> viewUp() const { return up; }
	const float viewTheta() const { return theta; }

private:
	CC_MODE			mode;

	Vector3f		eye;
	Vector3f		ref;
	Vector3f		up;
	float			theta;

	Vector3f		lastPointOnHemisphere;
	Vector3f		lastEye;
	Vector3f		lastRef;
	Vector3f		lastUp;
	int				lastTransX;
	int				lastTransY;
	int				lastPanY;

	int				screenWidth;
	int				screenHeight;

	Vector3f		eyeInitial;
	Vector3f		refInitial;
	Vector3f		upInitial;
	float			thetaInitial;
};

};	// namespace virtualspace

#endif