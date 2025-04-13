#ifndef _ILLUMINATION_SCENE_H_
#define _ILLUMINATION_SCENE_H_

#include "../../../Basis/Basis/Array.h"
#include "../../../Basis/Timer/Timer.h"
#include "../Core/Camera.h"
#include "../Core/Sampler.h"
#include "../Core/Film.h"
#include "../Primitive/Aggregate.h"
#include "../Light/Light.h"
#include "../Integrator/SurfaceIntegrator.h"
#include "../Integrator/VolumeIntegrator.h"

namespace illuminationspace {

using basis::Array;

class ILLUMINATIONDLL Scene	{
public : 
	//////////////////////////////////////////////////////////////////////////
	// Constructor/Destructor
	Scene() { init(); }
	~Scene();

	void init();

	//////////////////////////////////////////////////////////////////////////
	// Setter	
	void setCamera(Camera *c)							{ camera = c;				}
	void setImageSampler(ImageSampler *sam)				{ imageSampler = sam;		}
	void setLights(std::vector<Light> lts)				{ lights.resize_set(lts);	}	
	void setBackground(Spectrum c)						{ background = c;			}
	void setAggregate(Aggregate agg)					{ aggregate = agg;			}
	void setIntegrator(Integrator *si)					{ integrator = si;			}	

#ifdef _USE_XML_
	//////////////////////////////////////////////////////////////////////////
	// XML IO
	static const std::string XmlElementName() { return std::string("Scene"); }
	void importXmlNode(xmlNodePtr node);	
#endif

	//////////////////////////////////////////////////////////////////////////
	// Bounding Box
	inline const BBox WorldBound() const;

	//////////////////////////////////////////////////////////////////////////
	// Rendering
	void render();
	void render_serial();
	void render_openmp();

	void integrate(Ray&, Spectrum&, const ImageSample&);

	inline const bool hit(const Ray& r, Intersection &intersection) const;
	inline const bool shadowHit(const Ray& r) const;

	inline const Spectrum Li(const Ray &r, real *alpha=NULL) const;

public :
	Camera*					camera;	
	ImageSampler*			imageSampler;

	Array<Light>			lights;
	Spectrum				background;
	
	BBox					bbox;
	Aggregate				aggregate;	

	Integrator*				integrator;
};

inline const bool Scene::hit(const Ray& r, Intersection &intersection) const {
	intersection.tHit = LARGE_i;
	return aggregate.hit(r, intersection);
}

inline const bool Scene::shadowHit(const Ray& r) const {
	return aggregate.shadowHit(r);
}

inline const BBox Scene::WorldBound() const {
	return aggregate.worldBound();
}

inline const Spectrum Scene::Li(const Ray &r, real *alpha) const {
	Spectrum Lo = integrator->Li(this, r, alpha);
	return Lo;
}

};		// illuminationspace

#endif