#ifndef _VIRTUALSPACE_MATERIAL_H_
#define _VIRTUALSPACE_MATERIAL_H_

#include "../../Principle.h"
#include "../Property/MatLight.h"

namespace virtualspace {

using basis::RGBA;

//////////////////////////////////////////////////////////////////////////
// Material Class Definition
//////////////////////////////////////////////////////////////////////////
class PRINCIPLEDLL Material : public MatLight {
public:
	//////////////////////////////////////////////////////////////////////////
	// Constructors	
	inline Material();
	inline Material(const Material &);
	inline Material(const RGBA &); 

	//////////////////////////////////////////////////////////////////////////
	// Assignment and set
	inline void init();
	inline const Material &set(const Material &);
	inline const Material &operator=(const Material &);	

	inline void setAllColor(const RGBA&);
	inline void setAlpha(const RGBA::value_type);	

	inline void setShininess(const RGBA::value_type s);
	inline void setEmissive(const RGBA &);
	inline void setEmissive(const RGBA::value_type r, const RGBA::value_type g, const RGBA::value_type b, const RGBA::value_type a);

	//////////////////////////////////////////////////////////////////////////
	// Getter
	inline const RGBA&				getEmissive()  const { return emissive;  }
	inline const RGBA::value_type	getShininess() const { return shininess; }

	//////////////////////////////////////////////////////////////////////////
	// Manipulation
	void glBindMaterial() const;
	void glBind_FrontMaterial() const;
	void glBind_BackMaterial() const;

	//////////////////////////////////////////////////////////////////////////
	// IO
	friend PRINCIPLEDLL std::ostream& operator<<(std::ostream& os, const Material& mat);

#ifdef _USE_XML_
	//////////////////////////////////////////////////////////////////////////
	// XML IO
	static inline const std::string XmlElementName() { return std::string("Material"); }
	void importXmlNode(xmlNodePtr node);
	const xmlNodePtr exportXmlNode() const;	
#endif

private:
	RGBA				emissive;
	RGBA::value_type	shininess;
};

//////////////////////////////////////////////////////////////////////////
// Constructors	
inline Material::Material() : MatLight(), shininess(0) {}
inline Material::Material(const Material& mat) : MatLight() {
	set(mat);
}

inline Material::Material(const RGBA& color) : MatLight() {
	ambient = diffuse = color;
}

//////////////////////////////////////////////////////////////////////////
// Assignment and set
inline void Material::init() {
	MatLight::init();
	emissive.init();
	shininess = 0.0f;
}

inline const Material &Material::set(const Material &mat) {
	ambient = mat.ambient;	
	diffuse = mat.diffuse;
	specular = mat.specular;	
	emissive = mat.emissive;
	shininess = mat.shininess;
	return (*this);
}

inline const Material &Material::operator=(const Material &mat) {
	return set(mat);
}

inline void Material::setAllColor(const RGBA& rgba) {
	emissive = rgba;
}

inline void Material::setAlpha(const RGBA::value_type alpha) {
	MatLight::setAlpha(alpha);
	emissive.setAlpha(alpha);
}

inline void Material::setShininess(const RGBA::value_type s) {	shininess = s;	}
inline void Material::setEmissive(const RGBA &c) { emissive = c; }
inline void Material::setEmissive(const RGBA::value_type r, const RGBA::value_type g, const RGBA::value_type b, const RGBA::value_type a) {
	emissive.set(r,g,b,a);
}

PRINCIPLEDLL std::ostream& operator<<(std::ostream& os, const Material& mat);

};		// virtualspace

#endif