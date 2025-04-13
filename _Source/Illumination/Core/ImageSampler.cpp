#include <stdafx.h>
#include <Basis/MemoryManagement/MemoryManagement.h>
#include <Illumination/IlluminationHeader/Core/ImageSampler.h>

using namespace illuminationspace;

//////////////////////////////////////////////////////////////////////////
// ImageSampler
ImageSampler::ImageSampler(int xStart, int xEnd, int yStart, int yEnd, int nSpp) {
	xPixelStart = xStart;
	xPixelEnd = xEnd;
	yPixelStart = yStart;
	yPixelEnd = yEnd;
	nSamplesPerPixel = nSpp;
}

#ifdef _USE_XML_

void ImageSampler::importXmlNode(xmlNodePtr node) {}

ImageSampler* ImageSampler::createFromXmlNode(xmlNodePtr node) {
	ImageSampler *ret = NULL;

	if(node == NULL)
		basis::throwError("No ImageSampler Node _ illumination::ImageSampler::createFromXmlNode");

	std::string type = basis::XmlLibrary::xmlGetString(node, "type");
	if(type.compare(StratifiedImageSampler::XmlElementName()) == 0)		
		ret = new StratifiedImageSampler(node);	

	return ret;
}

#endif

//////////////////////////////////////////////////////////////////////////
// Stratified
StratifiedImageSampler::StratifiedImageSampler(int xStart, int xEnd, int yStart, int yEnd, int nXs, int nYs, bool jitter)
: ImageSampler(xStart, xEnd, yStart, yEnd, nXs * nYs), imageSamples(NULL), lensSamples(NULL), timeSamples(NULL) {
	jitterSamples = jitter;	
	nXPixelSamples = nXs;
	nYPixelSamples = nYs;
	
	init(xStart, xEnd, yStart, yEnd);	
}

void StratifiedImageSampler::release() {
	if(imageSamples) 
		basis::FreeAligned(imageSamples);
	imageSamples = lensSamples = timeSamples = NULL;
}

void StratifiedImageSampler::init(int xStart, int xEnd, int yStart, int yEnd) {
	release();

	xPixelStart = xStart;
	xPixelEnd = xEnd;
	yPixelStart = yStart;
	yPixelEnd = yEnd;

	xPos = xPixelStart;
	yPos = yPixelStart;	

	imageSamples = (real*)basis::AllocAligned(5 * nXPixelSamples * nYPixelSamples * sizeof(real));
	lensSamples = imageSamples + 2 * nXPixelSamples * nYPixelSamples;
	timeSamples = lensSamples + 2 * nXPixelSamples * nYPixelSamples;

	SamplerLibrary::StratifiedSample2D(imageSamples, nXPixelSamples, nYPixelSamples, jitterSamples);
	SamplerLibrary::StratifiedSample2D(lensSamples, nXPixelSamples, nYPixelSamples, jitterSamples);
	SamplerLibrary::StratifiedSample1D(timeSamples, nXPixelSamples * nYPixelSamples, jitterSamples);

	for(int i=0;i<2*nXPixelSamples*nYPixelSamples;i+=2) {
		imageSamples[i] += xPos;
		imageSamples[i+1] += yPos;
	}

	SamplerLibrary::Shuffle(lensSamples, nXPixelSamples * nYPixelSamples, 2);
	SamplerLibrary::Shuffle(timeSamples, nXPixelSamples * nYPixelSamples, 1);

	idxSamplePos = 0;
}

const bool StratifiedImageSampler::getNextSample(ImageSample *sample) {
	if(idxSamplePos == nXPixelSamples * nYPixelSamples) {
		if(++xPos == xPixelEnd) {
			xPos = xPixelStart;
			++yPos;
		}
		if(yPos == yPixelEnd) {
			//xPos = xPixelStart;	yPos = yPixelStart;
			return false;
		}

		SamplerLibrary::StratifiedSample2D(imageSamples, nXPixelSamples, nYPixelSamples, jitterSamples);
		SamplerLibrary::StratifiedSample2D(lensSamples, nXPixelSamples, nYPixelSamples, jitterSamples);
		SamplerLibrary::StratifiedSample1D(timeSamples, nXPixelSamples * nYPixelSamples, jitterSamples);

		for(int i=0;i<2*nXPixelSamples*nYPixelSamples;i+=2) {
			imageSamples[i] += xPos;
			imageSamples[i+1] += yPos;
		}

		SamplerLibrary::Shuffle(lensSamples, nXPixelSamples * nYPixelSamples, 2);
		SamplerLibrary::Shuffle(timeSamples, nXPixelSamples * nYPixelSamples, 1);

		idxSamplePos = 0;   
	}

	sample->imageX = imageSamples[2*idxSamplePos];
	sample->imageY = imageSamples[2*idxSamplePos+1];
	sample->lensU = lensSamples[2*idxSamplePos];
	sample->lensV = lensSamples[2*idxSamplePos+1];
	sample->time = timeSamples[idxSamplePos];

	++idxSamplePos;

	return true;
}

void StratifiedImageSampler::getSamples(Array<ImageSample>& samples, int x, int y) const {

	const real dx = 1.0f/static_cast<real>(nXPixelSamples), dy = 1.0f/static_cast<real>(nYPixelSamples);
	
	Array<ImageSample>::index_type idx = 0;

	for(int i=0;i<nXPixelSamples;++i) {
		for(int j=0;j<nYPixelSamples;++j) {
			const real jx = jitterSamples ? random() : 0.5f;
			const real jy = jitterSamples ? random() : 0.5f;
			samples[idx  ].imageX = (i+jx) * dx;
			samples[idx++].imageY = (j+jy) * dy;
		}
	}

	for(size_t i=0;i<samples.size();++i) {
		samples[i].imageX += x;
		samples[i].imageY += y;
	}
}

#ifdef _USE_XML_

void StratifiedImageSampler::importXmlNode(xmlNodePtr node) {
	xmlNodePtr nodePixel = basis::XmlLibrary::xmlFindNode(node, "NumPixelSamples");

	nXPixelSamples = basis::XmlLibrary::xmlGetInt(nodePixel, "x");
	nYPixelSamples = basis::XmlLibrary::xmlGetInt(nodePixel, "y");

	jitterSamples = basis::XmlLibrary::xmlGetBool(node, "Jitter");

	nSamplesPerPixel = nXPixelSamples * nYPixelSamples;
}

#endif