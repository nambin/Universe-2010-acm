#include <stdafx.h>
#include <Principle/VIrtualSpace/Camera/Camera.h>
#include <iostream>

using namespace std;
using namespace arithmetic;
using namespace virtualspace;

Camera::Camera() {
	set(Vector3f(0,0,2), Vector3f(0,0,0), Vector3f(0,1,0), 45.0f);
	ioHandler = new Camera_IOHandler(this);
}

Camera::Camera(Vector3f _eye, Vector3f _ref, Vector3f _up, float _theta) {
	set(_eye, _ref, _up, _theta);
	ioHandler = new Camera_IOHandler(this);
}

void Camera::set() {
	lastEye	= eye = eyeInitial;
	lastRef	= ref = refInitial;
	lastUp	= up  = upInitial;
	theta	= thetaInitial;

	mode = CC_IDLE_MODE;	
}

void Camera::set(Vector3f _eye, Vector3f _ref, Vector3f _up, float _theta) {
	eyeInitial = lastEye = eye = _eye;
	refInitial = lastRef = ref = _ref;
	upInitial  = lastUp  = up  = _up;
	thetaInitial = theta = _theta;

	perspective();
	viewport();

	mode = CC_IDLE_MODE;
}

void Camera::glBind() const {
	viewport();

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	lookAt();
}

void Camera::lookAt() const {
	gluLookAt(eye[0], eye[1], eye[2], ref[0], ref[1], ref[2], up[0], up[1], up[2]);
}

void Camera::resizeScreen(int _width, int _height) {
	screenWidth		= _width;
	screenHeight	= _height;
	perspective();
	viewport();
}

void Camera::viewport() const {
	glViewport(0, 0, (GLsizei)screenWidth, (GLsizei)screenHeight);
}

const Matrix4f Camera::getModelViewMatrix() const {
	double mat[16];

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	lookAt();
	glGetDoublev(GL_MODELVIEW_MATRIX, mat);

	glPopMatrix();
 
	Matrix4f cmat;
	for(Matrix4f::index_type i=0;i<16;++i) cmat(i) = (float)mat[i];

	return cmat;
}

void Camera::perspective() const {
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(theta, (GLdouble)screenWidth/(GLdouble)screenHeight, PerspectiveNear(), PerspectiveFar());
}

const Camera::Ray Camera::getLastRay(int x, int y) const {
	Vector3f	clickedPt, eyeToRef, dir, ortho;
	float		len;

	eyeToRef	= lastRef - lastEye;
	ortho		= cross(eyeToRef, lastUp);
	ortho.normalize();

	len = eyeToRef.mag() * (float)tan(theta*0.5*DegToRad);

	clickedPt	=	lastRef;
	clickedPt	+=	lastUp	* len *	(1.0f - (float)y/((float)screenHeight*0.5f));
	clickedPt	+=	ortho	* len *	((float)x/((float)screenWidth*0.5f) - 1.0f) * ((float)screenWidth/(float)screenHeight);

	dir = clickedPt - lastEye;
	dir.normalize();

	return Camera::Ray(lastEye, dir);
}

const Camera::Ray Camera::getRay(int x, int y) const {
	Vector3f	clickedPt, eyeToRef, dir, ortho;
	float		len;

	eyeToRef	= ref - eye;
	ortho		= cross(eyeToRef, up);
	ortho.normalize();

	len = eyeToRef.mag() * (float)tan(theta*0.5*DegToRad);

	clickedPt	=	ref;
	clickedPt	+=	up		* len *	(1.0f - (float)y/((float)screenHeight*0.5f))	;
	clickedPt	+=	ortho	* len *	((float)x/((float)screenWidth*0.5f) - 1.0f) * ((float)screenWidth/(float)screenHeight);

	dir = clickedPt - eye;
	dir.normalize();

	return Camera::Ray(eye, dir);
}

const Vector3f Camera::hemisphereMap(int x, int y) const {
	Vector3f	org, dir, rtn;
	float		t, a, b, c;
	Camera::Ray	ray;

	dir = getLastRay(x, y).dir;
	org	= unitVector(lastEye-lastRef)*(1.0f/sin(theta*0.5f*(float)DegToRad));

	a	= dot(dir, dir);
	b	= dot(org, dir);
	c	= dot(org, org) - 1;

	if(SQ(b)-a*c<=0) {
		t	= -b;
	}
	else {
		t	= -b-sqrt(SQ(b)-a*c);
	}

	return  unitVector(org+t*dir);
}

void Camera::trackballMouseDn(int x, int y){
	lastEye						= eye;
	lastRef						= ref;
	lastUp						= up;
	lastPointOnHemisphere		= hemisphereMap(x, y);

	mode = CC_TRACKBALL_MODE;
}

void Camera::trackballMouseUp(int x, int y){
	mode = CC_IDLE_MODE;
}

void Camera::trackballMouseMove(int x, int y){
	Quaternion<float>	eyeToRefQ(lastEye - lastRef);		// Quaternion used as Vector3
	Quaternion<float>	upQ(lastUp);						// Quaternion used as Vector3
	Vector3f			currentPointOnHemiphere = hemisphereMap(x, y);
	Vector3f			axis = unitVector(cross(lastPointOnHemisphere, currentPointOnHemiphere));
	float				theta = (lastPointOnHemisphere - currentPointOnHemiphere).mag()*RotationSpeed();

	Quaternion<float> q(axis, -theta);

	eyeToRefQ	= q * eyeToRefQ * inverse(q);
	eye			= lastRef + eyeToRefQ.vec();
	up			= (q * upQ * inverse(q)).vec();
}

void Camera::translateMouseDn(int x, int y){
	lastEye		= eye;
	lastRef		= ref;
	lastUp		= up;
	lastTransX	= x;
	lastTransY	= y;

	mode		= CC_TRANS_MODE;
}

void Camera::translateMouseUp(int x, int y){
	lastEye		= eye;
	lastRef		= ref;
	lastUp		= up;

	mode		= CC_IDLE_MODE;
}

void Camera::translateMouseMove(int x, int y) {
	Vector3f othor	= cross(unitVector(lastEye-lastRef), lastUp);
	Vector3f mov	= (float)(x-lastTransX)*othor + (float)(y-lastTransY)*lastUp;
	mov		*=	(lastEye-lastRef).mag() * TransSpeed();
	ref = lastRef + mov;
	eye = lastEye + mov;
}

void Camera::panMouseDn(int x, int y) {
	lastEye		= eye;
	lastRef		= ref;
	lastUp		= up;
	lastPanY	= y;
	mode		= CC_PAN_MODE;
}

void Camera::panMouseUp(int x, int y) {
	lastEye		= eye;
	lastRef		= ref;
	lastUp		= up;
	mode		= CC_IDLE_MODE;
}

void Camera::panMouseMove(int x, int y) {
	pan(0.03f * PanSpeed() * (y-lastPanY));
	lastPanY = y;
}

void Camera::mouseMove(int x, int y) {
	switch(mode) {
		case CC_TRACKBALL_MODE :	trackballMouseMove(x, y);	break;
		case CC_TRANS_MODE :		translateMouseMove(x, y);	break;
		case CC_PAN_MODE :			panMouseMove(x, y);			break;
		case CC_IDLE_MODE :
		default:
			break;
	}
}

void Camera::pan(float mag) {
	Vector3f moving = (eye-ref) * mag;
	eye += moving;
//	cout << eye << endl;
}

void Camera::panin()	{ pan(-PanSpeed()); }
void Camera::panout()	{ pan(+PanSpeed()); }	

void Camera::printCamera() const {
	cout << "Eye : " << eye << endl;
	cout << "Ref : " << ref << endl;
	cout << "Up  : " << up << endl;
	cout << "Thetat : " << theta << endl;
}