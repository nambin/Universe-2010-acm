#ifndef _ILLUMINATION_IMAGE_SAMPLER_H_
#define _ILLUMINATION_IMAGE_SAMPLER_H_

#include <vector>
#include "../../../Basis/XmlLibrary/XmlLibrary.h"
#include "../Illumination/Random.h"
#include "../../IlluminationSpace.h"
#include "../Illumination/Constants.h"
#include "Sampler.h"

namespace illuminationspace {

using basis::RandomGen;

class Scene;
class Integrator;

//////////////////////////////////////////////////////////////////////////
// Image Sampling
struct ImageSample {
public : 
	real imageX, imageY, lensU, lensV, time;
};

class ILLUMINATIONDLL ImageSampler	{
public :
	ImageSampler() {}
	ImageSampler(int xStart, int xEnd, int yStart, int yEnd, int nSpp);

	virtual void init(int xStart, int xEnd, int yStart, int yEnd) = 0;
	virtual const bool getNextSample(ImageSample *s) = 0;
	virtual void getSamples(Array<ImageSample>&, int x, int y) const = 0;

	const int totalSamples() const {
		return nSamplesPerPixel * (xPixelEnd-xPixelStart) * (yPixelEnd-yPixelStart);
	}

#ifdef _USE_XML_
	//////////////////////////////////////////////////////////////////////////
	// XML IO
	static const std::string XmlElementName() { return std::string("ImageSampler"); }
	virtual void importXmlNode(xmlNodePtr node);

	static ImageSampler* createFromXmlNode(xmlNodePtr node);
#endif

public :
	int xPixelStart, xPixelEnd, yPixelStart, yPixelEnd;
	int nSamplesPerPixel;
};

class ILLUMINATIONDLL StratifiedImageSampler : public ImageSampler {
public :
	StratifiedImageSampler() { imageSamples = lensSamples = timeSamples = NULL; }	
	StratifiedImageSampler(int xStart, int xEnd, int yStart, int yEnd, int nXs, int nYs, bool jitter);

#ifdef _USE_XML_
	StratifiedImageSampler(xmlNodePtr node) { imageSamples = lensSamples = timeSamples = NULL; importXmlNode(node); }
#endif

	void release();

	void init(int xStart, int xEnd, int yStart, int yEnd);
	const bool getNextSample(ImageSample *s);	

	void getSamples(Array<ImageSample>&, int x, int y) const;

#ifdef _USE_XML_
	//////////////////////////////////////////////////////////////////////////
	// XML IO
	static const std::string XmlElementName() { return std::string("StratifiedImageSampler"); }
	virtual void importXmlNode(xmlNodePtr node);
#endif

private :
	int		nXPixelSamples, nYPixelSamples;
	bool	jitterSamples;

	real*	imageSamples, *lensSamples, *timeSamples;	
	int		idxSamplePos;
	int		xPos, yPos;
};

};		// illuminationspace

#endif