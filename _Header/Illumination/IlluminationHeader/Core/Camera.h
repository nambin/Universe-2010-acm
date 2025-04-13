#ifndef _ILLUMINATION_CAMERA_H_
#define _ILLUMINATION_CAMERA_H_

#include "../../../Principle/Arithmetic/Vector/Vector3.h"
#include "Ray.h"
#include "Film.h"
#include "Transform.h"
#include "Sampler.h"

namespace illuminationspace {

using arithmetic::Vector3;

class ILLUMINATIONDLL Camera {
public :	
	Camera() : film(NULL) {}
	Camera(const Transform &w2c, real hither, real yon, real sopen, real sclose, Film *f);
	~Camera();

	virtual const real generateRay(const ImageSample &sample, Ray *ray) const = 0;
	virtual void getPixel(real &a, real &b, const Vector3<real> &pt) const = 0;
	virtual const real projectedPixelSize(const Vector3<real> &pt) const = 0;

#ifdef _USE_XML_
	//////////////////////////////////////////////////////////////////////////
	// XML IO
	static const std::string XmlElementName() { return std::string("Camera"); }
	virtual void importXmlNode(xmlNodePtr node);

	static Camera* createFromXmlNode(xmlNodePtr node);
#endif

public :
	Transform	WorldToCamera, CameraToWorld;
	real		clipHither, clipYon;
	real		shutterOpen, shutterClose;

	Film *		film;
};

class ILLUMINATIONDLL ProjectiveCamera : public Camera {
public : 	
	ProjectiveCamera() {}
	ProjectiveCamera(const Transform &w2c, const Transform &proj, const real screen[4], real hither, real yon, real sopen, real sclose, real lensr, real focald, Film *f);	

	void getPixel(real &a, real &b, const Vector3<real> &pt) const;
	const real projectedPixelSize(const Vector3<real> &pt) const;

#ifdef _USE_XML_
	//////////////////////////////////////////////////////////////////////////
	// XML IO	
	static const std::string XmlElementName() { return std::string("ProjectiveCamera"); }
	void importXmlNode(xmlNodePtr node);
#endif

public :
	Transform	CameraToScreen, WorldToScreen, RasterToCamera;
	Transform	ScreenToRaster, RasterToScreen;
	real		lensRadius, focalDistance;

	real		projectedPixelSizeAtNear, projectedPixelSizeAtFar;
	real		zNearCamera, zFarCamera;

private :
	void setProjectedPixelSize();
};

class ILLUMINATIONDLL OrthoCamera : public ProjectiveCamera {
public : 		
	OrthoCamera(const Transform &w2c, const real screen[4], real hither, real yon, real sopen, real sclose, real lensr, real focald, Film *f);
	const real generateRay(const ImageSample &sample, Ray *ray) const;	

#ifdef _USE_XML_
	OrthoCamera(xmlNodePtr node);

	//////////////////////////////////////////////////////////////////////////
	// XML IO
	static const std::string XmlElementName() { return std::string("OrthoCamera"); }
	virtual void importXmlNode(xmlNodePtr node);
#endif

};

class ILLUMINATIONDLL PerspectiveCamera : public ProjectiveCamera {
public : 		
	PerspectiveCamera(const Transform &w2c, const real screen[4], real hither, real yon, real sopen, real sclose, real lensr, real focald, real fov, Film *f);
	const real generateRay(const ImageSample &sample, Ray *ray) const;	

#ifdef _USE_XML_
	PerspectiveCamera(xmlNodePtr node);

	//////////////////////////////////////////////////////////////////////////
	// XML IO
	static const std::string XmlElementName() { return std::string("PerspectiveCamera"); }
	void importXmlNode(xmlNodePtr node);
#endif
};

};		// illuminationspace

#endif