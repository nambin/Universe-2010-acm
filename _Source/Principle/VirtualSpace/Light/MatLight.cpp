#include <stdafx.h>
#include <Principle/VirtualSpace/Property/MatLight.h>

using namespace virtualspace;

#ifdef _USE_XML_

void MatLight::importXmlNode(xmlNodePtr node) {

	xmlNodePtr childNode = NULL;

	if(childNode = basis::XmlLibrary::xmlFindNode(node, "Ambient"))	ambient.importXmlNode(childNode);
	else															ambient.black();

	if(childNode = basis::XmlLibrary::xmlFindNode(node, "Diffuse"))	diffuse.importXmlNode(childNode);
	else															diffuse.black();

	if(childNode = basis::XmlLibrary::xmlFindNode(node, "Specular"))specular.importXmlNode(childNode);
	else															specular.black();
}

const xmlNodePtr MatLight::exportXmlNode() const {

	xmlNodePtr node = xmlNewNode(NULL, BAD_CAST "MatLight");

	if(!ambient.isBlack())	basis::XmlLibrary::xmlAddChild(node, "Ambient", ambient);
	if(!diffuse.isBlack())	basis::XmlLibrary::xmlAddChild(node, "Diffuse", diffuse);
	if(!specular.isBlack())	basis::XmlLibrary::xmlAddChild(node, "Specular", specular);

	return node;
}

#endif