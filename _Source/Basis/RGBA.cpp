#include <stdafx.h>
#include <Basis/Color/RGBA.h>

using namespace std;

namespace basis {

void RGBA::clamp(RGBA::value_type min, RGBA::value_type max) {
	for(size_t i=0;i<3;++i)	{
		if(data[i] > max)	data[i] = max;
		if(data[i] < min)	data[i] = min;
	}
}

#ifdef _USE_XML_

void RGBA::importXmlNode(xmlNodePtr node) {	
	if(node == NULL)
		init();

	node = basis::XmlLibrary::xmlFindNode(node, RGBA::XmlElementName().c_str());

	if(node == NULL)
		init();

	data[0] = basis::XmlLibrary::xmlGetFloat(node, "R");
	data[1] = basis::XmlLibrary::xmlGetFloat(node, "G");
	data[2] = basis::XmlLibrary::xmlGetFloat(node, "B");
	data[3] = basis::XmlLibrary::xmlGetFloat(node, "A");
}

const xmlNodePtr RGBA::exportXmlNode() const {

	xmlNodePtr node = xmlNewNode(NULL, BAD_CAST RGBA::XmlElementName().c_str());

	basis::XmlLibrary::xmlPutProperty(node, "R", r());
	basis::XmlLibrary::xmlPutProperty(node, "G", g());
	basis::XmlLibrary::xmlPutProperty(node, "B", b());
	basis::XmlLibrary::xmlPutProperty(node, "A", a());

	return node;
}

#endif

};		// virtualspace