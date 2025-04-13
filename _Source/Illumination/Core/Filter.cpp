#include <stdafx.h>
#include <Illumination/IlluminationHeader/Core/Filter.h>

using namespace illuminationspace;

#ifdef _USE_XML_

void Filter::importXmlNode(xmlNodePtr node) {
	xmlNodePtr widthNode = basis::XmlLibrary::xmlFindNode(node, "Width");
	if(widthNode == NULL)
		basis::throwError("No Width in Filter _ Filter::importXmlNode");

	xWidth = basis::XmlLibrary::xmlGetFloat(widthNode, "x");
	yWidth = basis::XmlLibrary::xmlGetFloat(widthNode, "y");
	invXWidth = 1.0f/xWidth;
	invYWidth = 1.0f/yWidth;
}

Filter* Filter::createFromXmlNode(xmlNodePtr node) {

	Filter *ret = NULL;

	if(node == NULL)
		basis::throwError("No Filter Node _ illumination::Filter::createFromXmlNode");

	std::string type = basis::XmlLibrary::xmlGetString(node, "type");
	if(type.compare(BoxFilter::XmlElementName()) == 0)				ret = new BoxFilter(node);
	else if(type.compare(TriangleFilter::XmlElementName()) == 0) 	ret = new TriangleFilter(node);		
	else if(type.compare(GaussianFilter::XmlElementName()) == 0)	ret = new GaussianFilter(node);
	else if(type.compare(MitchellFilter::XmlElementName()) == 0)	ret = new MitchellFilter(node);
	else if(type.compare(LanczosSincFilter::XmlElementName()) == 0)	ret = new LanczosSincFilter(node);
	else
		basis::throwError("No Right Filter Name _ Filter::createFromXmlNode");

	return ret;

}

void BoxFilter::importXmlNode(xmlNodePtr node) {
	Filter::importXmlNode(node);
}

void TriangleFilter::importXmlNode(xmlNodePtr node) {
	Filter::importXmlNode(node);
}

void GaussianFilter::importXmlNode(xmlNodePtr node) {
	Filter::importXmlNode(node);

	alpha = basis::XmlLibrary::xmlGetFloat(node, "Alpha");
	xmlNodePtr expNode = basis::XmlLibrary::xmlFindNode(node, "Exp");
	if(expNode == NULL)
		basis::throwError("No Exp in GaussianFilter _ GaussianFilter::importXmlNode");

	expX = basis::XmlLibrary::xmlGetFloat(expNode, "x");
	expY = basis::XmlLibrary::xmlGetFloat(expNode, "y");
}

void MitchellFilter::importXmlNode(xmlNodePtr node) {
	Filter::importXmlNode(node);

	B = basis::XmlLibrary::xmlGetFloat(node, "B");
	C = basis::XmlLibrary::xmlGetFloat(node, "C");
}

void LanczosSincFilter::importXmlNode(xmlNodePtr node) {
	Filter::importXmlNode(node);

	tau = basis::XmlLibrary::xmlGetFloat(node, "Tau");
}

#endif