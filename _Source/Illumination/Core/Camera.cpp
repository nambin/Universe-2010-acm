#include <stdafx.h>
#include <Principle/Arithmetic/MathLibrary/Interpolation.h>
#include <Illumination/IlluminationHeader/Core/Camera.h>
#include <Illumination/IlluminationHeader/Core/IlluminationLibrary.h>

using namespace illuminationspace;

//////////////////////////////////////////////////////////////////////////
// Camera
Camera::Camera(const Transform &w2c, real hither, real yon, real sopen, real sclose, Film *f) 
: clipHither(hither), clipYon(yon), shutterOpen(sopen), shutterClose(sclose), film(f) {
	WorldToCamera.set(w2c);
	CameraToWorld.set(w2c.getInverse());
}

Camera::~Camera() {
	if(film) delete film;
}

#ifdef _USE_XML_

Camera* Camera::createFromXmlNode(xmlNodePtr node) {

	Camera *ret = NULL;

	if(node == NULL)
		throwError("No Camera Node _ illumination::Camera::createFromXmlNode");

	std::string type = basis::XmlLibrary::xmlGetString(node, "type");
	if(type.compare(OrthoCamera::XmlElementName()) == 0)		
		ret = new OrthoCamera(node);
	else if(type.compare(PerspectiveCamera::XmlElementName()) == 0) 
		ret = new PerspectiveCamera(node);		

	return ret;
}

void Camera::importXmlNode(xmlNodePtr node) {

	xmlNodePtr w2cNode = basis::XmlLibrary::xmlFindNode(node, "World2Camera");
	if(w2cNode != NULL) {
		WorldToCamera.importXmlNode(w2cNode);
		CameraToWorld.set(WorldToCamera.getMatrixInvRef(), WorldToCamera.getMatrixRef());
	}
	else {
		WorldToCamera.identity();
		CameraToWorld.identity();
	}

	clipHither = static_cast<real>(basis::XmlLibrary::xmlGetFloat(node, "Hither"));
	clipYon = static_cast<real>(basis::XmlLibrary::xmlGetFloat(node, "Yon"));
	shutterOpen = static_cast<real>(basis::XmlLibrary::xmlGetFloat(node, "ShutterOpen"));
	shutterClose = static_cast<real>(basis::XmlLibrary::xmlGetFloat(node, "ShutterClose"));

	xmlNodePtr filmNode = basis::XmlLibrary::xmlFindNode(node, "Film");
	film = Film::createFromXmlNode(filmNode);
}

#endif

//////////////////////////////////////////////////////////////////////////
// Projective Camera
ProjectiveCamera::ProjectiveCamera(const Transform &w2c, const Transform &proj, const real screen[4], real hither, real yon, real sopen, real sclose, real lensr, real focald, Film *f)
: Camera(w2c, hither, yon, sopen, sclose, f) {

	lensRadius = lensr; 
	focalDistance = focald;
	CameraToScreen = proj;		

	WorldToScreen = CameraToScreen * WorldToCamera;
	ScreenToRaster =  arithmetic::scale(real(film->xResolution), real(film->yResolution), real(1.0f)) 
					* arithmetic::scale(real(1.0f) / (screen[1]-screen[0]), real(1.0f) / (screen[3] - screen[2]), real(1.0f)) 
					* arithmetic::translate(Vector3<real>(-screen[0], -screen[2], real(0.0f)));

	RasterToScreen.set(ScreenToRaster.getInverse());
	RasterToCamera = CameraToScreen.getInverse() * RasterToScreen;

	setProjectedPixelSize();
}

void ProjectiveCamera::setProjectedPixelSize() {
	Vector3<real> v0(0.0f, 0.0f, 0.0f), vx(1.0f, 0.0f, 0.0f), vy(0.0f, 1.0f, 0.0f);
	Vector3<real> v0_C, vx_C, vy_C;
	real pixelSizeX, pixelSizeY;

	v0.z() = vx.z() = vy.z() = -1.0f;
	v0_C = RasterToCamera.transformPoint(v0);
	vx_C = RasterToCamera.transformPoint(vx);
	vy_C = RasterToCamera.transformPoint(vy);

	pixelSizeX = length(v0_C, vx_C);
	pixelSizeY = length(v0_C, vy_C);

	projectedPixelSizeAtNear = MIN(pixelSizeX, pixelSizeY);

	v0.z() = vx.z() = vy.z() = 1.0f;
	v0_C = RasterToCamera.transformPoint(v0);
	vx_C = RasterToCamera.transformPoint(vx);
	vy_C = RasterToCamera.transformPoint(vy);

	pixelSizeX = length(v0_C, vx_C);
	pixelSizeY = length(v0_C, vy_C);

	projectedPixelSizeAtFar = MIN(pixelSizeX, pixelSizeY);
}

void ProjectiveCamera::getPixel(real &a, real &b, const Vector3<real> &pt) const {
	Vector3<real> pScreen = WorldToScreen.transformPoint(pt);
	Vector3<real> pRaster = ScreenToRaster.transformPoint(pScreen);
	a = pRaster[0];
	b = pRaster[1];
}

const real ProjectiveCamera::projectedPixelSize(const Vector3<real> &pt) const {
	Vector3<real> pCamera = WorldToCamera.transformPoint(pt);		

	return (  projectedPixelSizeAtNear * (pCamera[2]  - zFarCamera) 
			+ projectedPixelSizeAtFar  * (zNearCamera - pCamera[2]) ) / (zNearCamera - zFarCamera);
}

#ifdef _USE_XML_

void ProjectiveCamera::importXmlNode(xmlNodePtr node) {

	Camera::importXmlNode(node);	

	node = basis::XmlLibrary::xmlFindNode(node, ProjectiveCamera::XmlElementName().c_str());
	if(node == NULL)
		throwError("No ProjectiveCamera Node _ illumination::ProjectiveCamera::importXmlNode");

	lensRadius = basis::XmlLibrary::xmlGetFloat(node, "LensRadius");
	focalDistance = basis::XmlLibrary::xmlGetFloat(node, "FocalDistance");

	real screen[4];
	xmlNodePtr screenNode = basis::XmlLibrary::xmlFindNode(node, "ScreenWindow");
	if(screenNode != NULL) {
		screen[0] = basis::XmlLibrary::xmlGetFloat(screenNode, "Left");
		screen[1] = basis::XmlLibrary::xmlGetFloat(screenNode, "Right");
		screen[2] = basis::XmlLibrary::xmlGetFloat(screenNode, "Down");
		screen[3] = basis::XmlLibrary::xmlGetFloat(screenNode, "Up");
	}
	else 
		throwError("No Screen Window _ illumination::ProjectiveCamera::importXmlNode");

	//////////////////////////////////////////////////////////////////////////
	// Modify Screen
	real vertical = screen[3] - screen[2], horizontal, mid_horizontal;
	real xResolution = static_cast<real>(film->xResolution), yResolution = static_cast<real>(film->yResolution);

	mid_horizontal = (screen[0] + screen[1]) * 0.5f;
	horizontal = vertical / yResolution * xResolution;
	screen[0] = mid_horizontal - horizontal * 0.5f;
	screen[1] = mid_horizontal + horizontal * 0.5f;

	//////////////////////////////////////////////////////////////////////////
	// Perspective Camera
	xmlNodePtr perspectiveNode = basis::XmlLibrary::xmlFindNode(node, PerspectiveCamera::XmlElementName().c_str());
	if(perspectiveNode != NULL) {
		real fov = basis::XmlLibrary::xmlGetFloat(perspectiveNode, "Fov");
		CameraToScreen = arithmetic::perspective(fov, clipHither, clipYon);
	}
	else {
		CameraToScreen = arithmetic::orthographic(clipHither, clipYon);
	}

	//////////////////////////////////////////////////////////////////////////
	// Setting
	WorldToScreen = CameraToScreen * WorldToCamera;
	ScreenToRaster =  arithmetic::scale(real(film->xResolution), real(film->yResolution), real(1.0f)) 
					* arithmetic::scale(real(1.0f) / (screen[1]-screen[0]), real(1.0f) / (screen[3] - screen[2]), real(1.0f)) 
					* arithmetic::translate(Vector3<real>(-screen[0], -screen[2], real(0.0f)));

	RasterToScreen.set(ScreenToRaster.getInverse());
	RasterToCamera = CameraToScreen.getInverse() * RasterToScreen;

	setProjectedPixelSize();
}

//////////////////////////////////////////////////////////////////////////
// Orthographic Camera
OrthoCamera::OrthoCamera(xmlNodePtr node) {
	importXmlNode(node);
}

#endif

OrthoCamera::OrthoCamera(const Transform &w2c, const real screen[4], real hither, real yon, real sopen, real sclose, real lensr, real focald, Film *f)
: ProjectiveCamera(w2c, arithmetic::orthographic(hither, yon), screen, hither, yon, sopen, sclose, lensr, focald, f) {
	zNearCamera = -clipHither;
	zFarCamera = -clipYon;
}

const real OrthoCamera::generateRay(const ImageSample &sample, Ray *ray) const {
	Vector3<real> pRas(sample.imageX, sample.imageY, 0.0f);
	Vector3<real> pCamera = RasterToCamera.transformPoint(pRas);

	ray->setOrigin(pCamera);
	ray->setDirection(Vector3<real>(0.0f, 0.0f, -1.0f));

	// set time, depth of field
	ray->time = arithmetic::lerp(shutterOpen, shutterClose, sample.time);
	if(lensRadius > 0.0f) {
		real lensU, lensV;
		IlluminationLibrary::ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
		lensU *= lensRadius;
		lensV *= lensRadius;
		real ft = (focalDistance - clipHither) / ray->direction().z();
		Vector3<real> pFocus = (*ray)(ft);

		(ray->origin())[0] += lensU;
		(ray->origin())[1] += lensV;
		ray->setDirection(pFocus - ray->origin());
	}

	ray->setTmin(0.0f);
	ray->setTmax(clipYon - clipHither);	
	CameraToWorld.transformRay(*ray, ray);
	//(*ray) = CameraToWorld.transformRay(*ray);

	return 1.0f;
}

#ifdef _USE_XML_

void OrthoCamera::importXmlNode(xmlNodePtr node) {	
	ProjectiveCamera::importXmlNode(node);
}

//////////////////////////////////////////////////////////////////////////
// Perspective Camera
PerspectiveCamera::PerspectiveCamera(xmlNodePtr node) {
	importXmlNode(node);
}

#endif

PerspectiveCamera::PerspectiveCamera(const Transform &w2c, const real screen[4], real hither, real yon, real sopen, real sclose, real lensr, real focald, real fov, Film *f)
: ProjectiveCamera(w2c, arithmetic::perspective(fov, hither, yon), screen, hither, yon, sopen, sclose, lensr, focald, f) {
	zNearCamera = -clipHither;
	zFarCamera = -clipYon;
}

const real PerspectiveCamera::generateRay(const ImageSample &sample, Ray *ray) const {
	Vector3<real> pRas(sample.imageX, sample.imageY, 0.0f);
	Vector3<real> pCamera = RasterToCamera.transformPoint(pRas);

	ray->setOrigin(pCamera);							//ray->setOrigin(Vector3<real>(0.0f));
	ray->setDirection(unitVector(pCamera));

	// set time, depth of field
	ray->time = arithmetic::lerp(shutterOpen, shutterClose, sample.time);
	if(lensRadius > 0.0f) {
		real lensU, lensV;
		IlluminationLibrary::ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
		lensU *= lensRadius;
		lensV *= lensRadius;
		real ft = (focalDistance - clipHither) / ray->direction().z();
		Vector3<real> pFocus = (*ray)(ft);

		(ray->origin())[0] += lensU;
		(ray->origin())[1] += lensV;
		ray->setDirection(pFocus - ray->origin());
	}

	ray->setTmin(0.0f);
	ray->setTmax((clipYon - clipHither) / ray->direction().z());	
	CameraToWorld.transformRay(*ray, ray);
	//(*ray) = CameraToWorld.transformRay(*ray);

	return 1.0f;
}

#ifdef _USE_XML_

void PerspectiveCamera::importXmlNode(xmlNodePtr node) {	
	ProjectiveCamera::importXmlNode(node);
}

#endif