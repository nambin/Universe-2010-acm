#include <stdafx.h>
#include <Nature/MassSpring/Core/ClothLibrary.h>

using namespace springspace;

const Material ClothLibrary::Material_Front_Basic() {
	RGBA c = RGBA(255, 15, 250) / 255;

	Material mat;
	mat.setAmbient(c);
	mat.setDiffuse(c);
	mat.setSpecular(c);
	mat.setEmissive(RGBA::black());
	mat.setShininess(0);
	return mat;
}

const Material ClothLibrary::Material_Back_Basic() {
	RGBA c = RGBA(160, 120, 80) / 255;

	Material mat;
	mat.setAmbient(c);
	mat.setDiffuse(c);
	mat.setSpecular(c);
	mat.setEmissive(RGBA::black());
	mat.setShininess(0);
	return mat;
}