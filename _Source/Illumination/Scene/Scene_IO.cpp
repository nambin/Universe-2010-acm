#include <stdafx.h>
#include <Basis/XmlLibrary/XmlLibrary.h>
#include <Illumination/IlluminationHeader/Scene/Scene.h>

using namespace illuminationspace;

#ifdef _USE_XML_
void Scene::importXmlNode(xmlNodePtr node) {	
	if(node == NULL)
		throwError("No Scene Node _ illumination::Scene::importXmlNode");

	xmlNodePtr sceneNode = basis::XmlLibrary::xmlFindNode(node, illuminationspace::Scene::XmlElementName().c_str());
	if(sceneNode == NULL)
		throwError("No Scene Node _ illumination::Scene::importXmlNode");

	xmlNodePtr backgroundNode = basis::XmlLibrary::xmlFindNode(node, "BackGroundColor");
	if(backgroundNode)	background.importXmlNode(backgroundNode);
	else				background.init();

	xmlNodePtr cameraNode = basis::XmlLibrary::xmlFindNode(sceneNode, illuminationspace::Camera::XmlElementName().c_str());
	camera = Camera::createFromXmlNode(cameraNode);

	xmlNodePtr imageSamplerNode = basis::XmlLibrary::xmlFindNode(sceneNode, ImageSampler::XmlElementName().c_str());
	imageSampler = ImageSampler::createFromXmlNode(imageSamplerNode);

	int xstart, xend, ystart, yend;
	camera->film->getSampleExtent(xstart, xend, ystart, yend);
	imageSampler->init(xstart, xend, ystart, yend);

}
#endif