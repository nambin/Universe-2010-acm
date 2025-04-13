#ifndef _VIRTUALSPACE_MATLIGHT_H_
#define _VIRTUALSPACE_MATLIGHT_H_

#include "../../Principle.h"
#include "../RenderingPipeline/RenderingPipelineSetting.h"
#include "../../../Basis/XmlLibrary/XmlLibrary.h"

namespace virtualspace {

using basis::RGBA;

class PRINCIPLEDLL MatLight	{
public:
	//////////////////////////////////////////////////////////////////////////
	// Constructors	
	inline MatLight();
	inline MatLight(const MatLight &);

	//////////////////////////////////////////////////////////////////////////
	// Assignment and set
	inline void init();
	inline void set(const MatLight &);

	inline void setAllColor(const RGBA&);
	inline void setAlpha(const RGBA::value_type);

	inline void setAmbient (const RGBA &);
	inline void setDiffuse (const RGBA &);
	inline void setSpecular(const RGBA &);	

	inline void setAmbient (RGBA::value_type r, RGBA::value_type g, RGBA::value_type b, RGBA::value_type a);
	inline void setDiffuse (RGBA::value_type r, RGBA::value_type g, RGBA::value_type b, RGBA::value_type a);
	inline void setSpecular(RGBA::value_type r, RGBA::value_type g, RGBA::value_type b, RGBA::value_type a);
	
	//////////////////////////////////////////////////////////////////////////
	// Getter 
	inline const RGBA &getAmbient () const { return ambient;  }
	inline const RGBA &getDiffuse () const { return diffuse;  }	
	inline const RGBA &getSpecular() const { return specular; }

#ifdef _USE_XML_
	//////////////////////////////////////////////////////////////////////////
	// XML IO
	void importXmlNode(xmlNodePtr node);
	const xmlNodePtr exportXmlNode() const;	
#endif

protected :
	RGBA	ambient,
			diffuse,
			specular;	
};

inline MatLight::MatLight() {}
inline MatLight::MatLight(const MatLight &v) {
	set(v);
}

inline void MatLight::init() {
	ambient.init();
	diffuse.init();
	specular.init();
}

inline void MatLight::set(const MatLight &v) {
	ambient = v.ambient;	
	diffuse = v.diffuse;
	specular = v.specular;	
}

inline void MatLight::setAllColor(const RGBA& rgba) {
	ambient = rgba;
	diffuse = rgba;
	specular = rgba;
}

inline void MatLight::setAlpha(const RGBA::value_type alpha) {
	ambient.setAlpha(alpha);
	diffuse.setAlpha(alpha);
	specular.setAlpha(alpha);
}

inline void MatLight::setAmbient(const RGBA &c) {	ambient = c;	}
inline void MatLight::setDiffuse(const RGBA &c) {	diffuse = c;	}
inline void MatLight::setSpecular(const RGBA &c) {	specular = c;	}

inline void MatLight::setAmbient(RGBA::value_type r, RGBA::value_type g, RGBA::value_type b, RGBA::value_type a) {	ambient.set(r,g,b,a);	}
inline void MatLight::setDiffuse(RGBA::value_type r, RGBA::value_type g, RGBA::value_type b, RGBA::value_type a) {	diffuse.set(r,g,b,a);	}
inline void MatLight::setSpecular(RGBA::value_type r, RGBA::value_type g, RGBA::value_type b, RGBA::value_type a){	specular.set(r,g,b,a);	}

};		// namespace virtualspace

#endif