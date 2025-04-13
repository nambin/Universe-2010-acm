#ifndef _ILLUMINATION_FILM_H_
#define _ILLUMINATION_FILM_H_

#include <cstdlib>
#ifdef _WIN32
#include "gl/glut.h"
#else
#include <IlluminationHeader/GLUT/glut.h>
#endif

#include <string>
#include "../../../Basis/Basis/Array.h"
#include "../../../Basis/Basis/Array2.h"
#include "Spectrum.h"
#include "Ray.h"
#include "Sampler.h"
#include "ImageSampler.h"
#include "Filter.h"
#include "../../IlluminationSpace.h"

namespace illuminationspace {

using std::string;
using basis::Array2;

class ILLUMINATIONDLL Film {
public : 
	Film() { /*nIter = 0;*/ }
	Film(int xres, int yres) : xResolution(xres), yResolution(yres) { /*nIter = 0;*/ }
	
	virtual void addSample(const ImageSample &, const Ray &, const Spectrum &, real alpha) = 0;
	virtual void writeImage(real gamma = 0.0f) = 0;
	virtual void getSampleExtent(int &xstart, int &xend, int &ystart, int &yend) const = 0;
	
	virtual void print_debug_info() const = 0;

#ifdef _USE_XML_
	//////////////////////////////////////////////////////////////////////////
	// XML IO
	static const std::string XmlElementName() { return std::string("Film"); }
	virtual void importXmlNode(xmlNodePtr node);

	static Film* createFromXmlNode(xmlNodePtr node);
#endif
	
public :
	int		xResolution, yResolution;
	real	gamma;

	//int nIter;
};

class ILLUMINATIONDLL ImageFilm : public Film {
public :
	struct Pixel {
		Pixel() : L(0.0f), alpha(0.0f), weightSum(0.0f) {};

		Spectrum L;
		real alpha, weightSum;
	};

public :
	//////////////////////////////////////////////////////////////////////////
	// ImageFilm Public Methods
	ImageFilm() : filter(NULL) {}	
	ImageFilm(int xres, int yres, Filter *filt, const real crop[4], const string &filename, bool premultiplyAlpha, int writingfrequency);
	~ImageFilm() {
		if(filter) delete filter;
	}

	void init();

	void addSample(const ImageSample &sample, const Ray &ray, const Spectrum &L, real alpha);
	void writeImage(real gamma=0.0f);
	void getSampleExtent(int &xstart, int &xend, int &ystart, int &yend) const;	

	void print_debug_info() const;

#ifdef _USE_XML_
	ImageFilm(xmlNodePtr node);

	//////////////////////////////////////////////////////////////////////////
	// XML IO
	static const std::string XmlElementName() { return std::string("ImageFilm"); }
	virtual void importXmlNode(xmlNodePtr node);
#endif

public : 	
	Filter*		filter;
	int			writeFrequency, sampleCount;
	string		filename;
	bool		premultiplyAlpha;
	real		cropWindow[4];
	int			xPixelStart, yPixelStart, xPixelCount, yPixelCount;

	Array2<Pixel>	pixels;
	Array2<real>	filterTable;
};

};		// illuminationspace

#endif
