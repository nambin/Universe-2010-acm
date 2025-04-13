#ifndef _ILLUMINATION_TEXTURE_H_
#define _ILLUMINATION_TEXTURE_H_

#include "../../../Basis/MemoryManagement/RCPtr.h"
#include "../../IlluminationSpace.h"
#include "../Shape/Shape.h"

namespace illuminationspace {

template <typename T> 
class Texture  {
protected :	
	class TextureValue;

public : 
	Texture() : value(NULL) {}
	Texture(TextureValue *ptr) : value(ptr) {}

	T evaluate(const DifferentialGeometry &dg) const { return value->evaluate(dg); }

protected :	
	class TextureValue : public basis::RCObject {
	public : 
		TextureValue() {}
		virtual ~TextureValue() {}

		TextureValue *clone() const = 0;

		virtual T evaluate(const DifferentialGeometry &) const = 0;
	};

	//////////////////////////////////////////////////////////////////////////
	// Reference Counted Pointer
	basis::RCPtr<TextureValue> value;

public : 
	const basis::RCPtr<TextureValue> &getTextureValue() const { return value; }	
};

template <typename T> 
class ConstantTexture : public Texture<T> {
public :
	ConstantTexture(const T& v) : Texture<T>(new ConstantTextureValue(v)) {}
	
protected :	
	class ConstantTextureValue : public Texture<T>::TextureValue {
	public : 
		ConstantTextureValue(const T &v) : value(v) {};

		ConstantTextureValue *clone() const { return new ConstantTextureValue(*this); }

		T evaluate(const DifferentialGeometry &) const {
			return value;
		}

		T value;
	};
};

};		// illuminationspace

#endif