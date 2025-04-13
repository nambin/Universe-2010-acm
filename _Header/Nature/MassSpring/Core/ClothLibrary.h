#ifndef _MASS_SPRING_CLOTH_LIBRARY_H_
#define _MASS_SPRING_CLOTH_LIBRARY_H_

#include "../../../Basis/Color/RGBA.h"
#include "../../../Principle/Virtualspace/Material/Material.h"

namespace springspace {

using basis::RGBA;
using virtualspace::Material;

class ClothLibrary {
public : 
	//////////////////////////////////////////////////////////////////////////
	// Material
	static const Material Material_Front_Basic();
	static const Material Material_Back_Basic();
};

};		// springspace

#endif