#include <stdafx.h>
#include <Illumination/IlluminationHeader/Core/Spectrum.h>

using namespace illuminationspace;

//////////////////////////////////////////////////////////////////////////
// Manipulation
const bool Spectrum::isBlack() const {
	for(int i=0;i<N_COLOR_SAMPLES;++i)
		if(c[i] != 0.0f)
			return false;
	return true;
}

void Spectrum::clamp(real min, real max)	{
	for(int i=0;i<N_COLOR_SAMPLES;++i)	{
		if(c[i] > max)	c[i] = max;
		if(c[i] < min)	c[i] = min;
	}
}

const Spectrum Spectrum::getClamped(real min, real max) const {
	Spectrum ret;
	for(int i=0;i<N_COLOR_SAMPLES;++i)	{
		if(c[i] < min)		ret.c[i] = min;
		else if(c[i] > max) ret.c[i] = max;
		else				ret.c[i] = c[i];
	}

	return ret;
}

const Spectrum& Spectrum::exponent()	{
	for(int i=0;i<N_COLOR_SAMPLES;++i)
		c[i] = ::exp(c[i]);
	return (*this);
}

const Spectrum& Spectrum::sqrt() {
	for(int i=0;i<N_COLOR_SAMPLES;++i)
		c[i] = std::sqrt(c[i]);
	return (*this);
}

#ifdef _USE_XML_
void Spectrum::importXmlNode(xmlNodePtr node) {
	if(node == NULL)
		init();

	node = basis::XmlLibrary::xmlFindNode(node, Spectrum::XmlElementName().c_str());

	if(node == NULL)
		init();

	c[0] = basis::XmlLibrary::xmlGetFloat(node, "R");
	c[1] = basis::XmlLibrary::xmlGetFloat(node, "G");
	c[2] = basis::XmlLibrary::xmlGetFloat(node, "B");
}
#endif