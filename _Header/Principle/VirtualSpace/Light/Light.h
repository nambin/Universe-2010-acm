#ifndef _VIRTUALSPACE_LIGHT_H_
#define _VIRTUALSPACE_LIGHT_H_

#include "../Property/MatLight.h"
#include "../../Arithmetic/Vector/Vector3.h"
#include "../../Arithmetic/Vector/Vector4.h"

namespace virtualspace {

using basis::RGBA;

//////////////////////////////////////////////////////////////////////////
// Light Class Definition
//////////////////////////////////////////////////////////////////////////
class PRINCIPLEDLL Light : public MatLight {
public:
	//////////////////////////////////////////////////////////////////////////
	// Static
	static int lightCount;
	static void initLights() { lightCount = 0; }

	//////////////////////////////////////////////////////////////////////////
	// Constructors	
	inline Light();
	inline Light(const Light &);

	//////////////////////////////////////////////////////////////////////////
	// Assignment and set
	inline const Light &set(const Light &);
	inline const Light &operator=(const Light &);	

	void     setAttenuation(const arithmetic::Vector3f &);
	void     setAttenuation(const float _const, const float _linear, const float _quad);

	void     setPosition(const arithmetic::Vector4f &);
	void     setPosition(const float x, const float y, const float z, const float w);

	void     setEnabled(bool);	

	//////////////////////////////////////////////////////////////////////////
	// Getter
	const arithmetic::Vector3f&	getAttenuation() const { return attenuation; }
	const arithmetic::Vector4f&	getPosition()    const { return position;    }
	const int					getLightID()     const { return lightID;     }

	const bool					isEnabled()		 const { return enabled;	 }

	//////////////////////////////////////////////////////////////////////////
	// Manipulation
	void	glDraw(float size = 1.1f) const;
	void	glBindLight() const;
	void	glUnbindLight() const;

private :
	arithmetic::Vector3f	attenuation;
	arithmetic::Vector4f	position; 
				
	bool		enabled;
	mutable int	lightID;
};

//////////////////////////////////////////////////////////////////////////
// Constructors	
inline Light::Light() : MatLight(), attenuation(arithmetic::Vector3f(1,0,0)), enabled(true), lightID(-1) {}
inline Light::Light(const Light &light) : MatLight() {
	set(light);
}

//////////////////////////////////////////////////////////////////////////
// Assignment and set
inline const Light &Light::set(const Light &light) {
	ambient = light.ambient;	
	diffuse = light.diffuse;
	specular = light.specular;	
	attenuation = light.attenuation;
	position = light.position;
	enabled = light.enabled;
	lightID = light.lightID;

	return (*this);
}

inline const Light &Light::operator=(const Light &light) {
	return set(light);
}

inline void Light::setAttenuation(const arithmetic::Vector3f &att) { attenuation = att; }
inline void Light::setAttenuation(float _const, float _linear, float _quad) {
	attenuation[0] = _const; attenuation[1] = _linear; attenuation[2] = _quad;
}

inline void Light::setPosition(const arithmetic::Vector4f &_pos) { position = _pos; }
inline void Light::setPosition(const float x, const float y, const float z, const float w) {
	position[0] = x; position[1] = y; position[2] = z; position[3] = w;
}

inline void Light::setEnabled(bool b)	{ enabled = b;	}

};		// virtualspace

#endif