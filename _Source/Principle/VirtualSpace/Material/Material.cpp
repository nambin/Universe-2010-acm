#include <stdafx.h>
#include <Principle/VirtualSpace/Material/Material.h>

using namespace virtualspace;

namespace virtualspace {

//////////////////////////////////////////////////////////////////////////
// Manipulation
void Material::glBindMaterial() const {
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT , ambient.getptr());
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE , diffuse.getptr());
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular.getptr());
	glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, emissive.getptr());
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
}

void Material::glBind_FrontMaterial() const {
	glMaterialfv(GL_FRONT, GL_AMBIENT , ambient.getptr());
	glMaterialfv(GL_FRONT, GL_DIFFUSE , diffuse.getptr());
	glMaterialfv(GL_FRONT, GL_SPECULAR, specular.getptr());
	glMaterialfv(GL_FRONT, GL_EMISSION, emissive.getptr());
	glMaterialf(GL_FRONT, GL_SHININESS, shininess);
}

void Material::glBind_BackMaterial() const {
	glMaterialfv(GL_BACK, GL_AMBIENT , ambient.getptr());
	glMaterialfv(GL_BACK, GL_DIFFUSE , diffuse.getptr());
	glMaterialfv(GL_BACK, GL_SPECULAR, specular.getptr());
	glMaterialfv(GL_BACK, GL_EMISSION, emissive.getptr());
	glMaterialf(GL_BACK, GL_SHININESS, shininess);
}

std::ostream& operator<<(std::ostream& os, const Material& mat) {
	os << "Ambient : " << mat.ambient << std::endl;
	os << "Diffuse : " << mat.diffuse << std::endl;
	os << "Specular : " << mat.specular << std::endl;
	os << "Emissive : " << mat.emissive << std::endl;
	os << "Shininess : " << mat.shininess << std::endl;
	return os;
}

#ifdef _USE_XML_

void Material::importXmlNode(xmlNodePtr node) {
	if(node == NULL)
		init();

	node = basis::XmlLibrary::xmlFindNode(node, Material::XmlElementName().c_str());
	if(node == NULL)
		init();

	MatLight::importXmlNode(node);

	xmlNodePtr childNode = NULL;

	if(childNode = basis::XmlLibrary::xmlFindNode(node, "Emissive"))	emissive.importXmlNode(childNode);
	else																emissive.black();

	if(childNode = basis::XmlLibrary::xmlFindNode(node, "Shininess"))	shininess = basis::XmlLibrary::xmlGetFloat(childNode, "Shininess");
	else																shininess = 0.0f;
}

const xmlNodePtr Material::exportXmlNode() const {
	xmlNodePtr node = MatLight::exportXmlNode();
	xmlNodeSetName(node, BAD_CAST "Material");

	if(!emissive.isBlack())	basis::XmlLibrary::xmlAddChild(node, "Emissive", emissive);
	if(shininess != 0.0f)	basis::XmlLibrary::xmlAddChild(node, "Shininess", shininess);

	return node;
}

#endif

};		// virtualspace