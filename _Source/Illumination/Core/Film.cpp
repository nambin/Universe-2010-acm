#include <stdafx.h>
#include <Basis/Image/Image.h>
#include <Illumination/IlluminationHeader/Core/Film.h>
#include <Illumination/IlluminationHeader/Illumination/BasicLibrary.h>

using namespace std;
using namespace illuminationspace;

#ifdef _USE_XML_

//////////////////////////////////////////////////////////////////////////
// Film
void Film::importXmlNode(xmlNodePtr node) {	
	xmlNodePtr resolutionNode = basis::XmlLibrary::xmlFindNode(node, "Resolution");

	xResolution = basis::XmlLibrary::xmlGetInt(resolutionNode, "x");
	yResolution = basis::XmlLibrary::xmlGetInt(resolutionNode, "y");
	gamma = basis::XmlLibrary::xmlGetFloat(node, "Gamma");
}

Film* Film::createFromXmlNode(xmlNodePtr node) {

	Film *ret = NULL;

	if(node == NULL)
		basis::throwError("No Film Node _ illumination::Film::createFromXmlNode");

	std::string type = basis::XmlLibrary::xmlGetString(node, "type");
	if(type.compare(ImageFilm::XmlElementName()) == 0)		
		ret = new ImageFilm(node);

	return ret;
}

//////////////////////////////////////////////////////////////////////////
// ImageFilm
ImageFilm::ImageFilm(xmlNodePtr node) : filter(NULL) {
	importXmlNode(node);
}

#endif

ImageFilm::ImageFilm(int xres, int yres, Filter *filt, const real crop[4], const string &fn, bool premult, int wf) 
: Film(xres, yres) {
	filter = filt;
	memcpy(cropWindow, crop, 4 * sizeof(real));
	filename = fn;
	premultiplyAlpha = premult;
	writeFrequency = sampleCount = wf;

	init();
}

void ImageFilm::init() {

#define FILTER_TABLE_SIZE 16

	xPixelStart = Ceil2Int(xResolution * cropWindow[0]);
	xPixelCount = Ceil2Int(xResolution * cropWindow[1]) - xPixelStart;
	yPixelStart = Ceil2Int(yResolution * cropWindow[2]);
	yPixelCount = Ceil2Int(yResolution * cropWindow[3]) - yPixelStart;

	pixels.resize(xPixelCount, yPixelCount);

	//////////////////////////////////////////////////////////////////////////
	// Pre_compute filter weight table
	filterTable.resize(FILTER_TABLE_SIZE, FILTER_TABLE_SIZE);
	
	for (int x = 0; x < FILTER_TABLE_SIZE; ++x) {
		real fx = ((real)x + 0.5f) * filter->xWidth / FILTER_TABLE_SIZE;
		for (int y = 0; y < FILTER_TABLE_SIZE; ++y) {
			real fy = ((real)y + 0.5f) * filter->yWidth / FILTER_TABLE_SIZE;
			filterTable(x,y) = filter->evaluate(fx, fy);
		}
	}
}

void ImageFilm::addSample(const ImageSample &sample, const Ray &ray, const Spectrum &L, real alpha) {

	//////////////////////////////////////////////////////////////////////////
	// Compute sample's raster extent
	real dImageX = sample.imageX - 0.5f;
	real dImageY = sample.imageY - 0.5f;
	int x0 = Ceil2Int(dImageX - filter->xWidth); int x1 = Floor2Int(dImageX + filter->xWidth);
	int y0 = Ceil2Int(dImageY - filter->yWidth); int y1 = Floor2Int(dImageY + filter->yWidth);

	x0 = MAX(x0, xPixelStart);	x1 = MIN(x1, xPixelStart + xPixelCount - 1);
	y0 = MAX(y0, yPixelStart);	y1 = MIN(y1, yPixelStart + yPixelCount - 1);

	//////////////////////////////////////////////////////////////////////////
	// Loop over filter support and add sample to pixel arrays
	// Pre_compute x and y filter table offsets
	int *ifx = (int *)alloca((x1-x0+1) * sizeof(int));
	for (int x = x0; x <= x1; ++x) {
		real fx = abs((x - dImageX) * filter->invXWidth * FILTER_TABLE_SIZE);
		ifx[x-x0] = MIN(Floor2Int(fx), FILTER_TABLE_SIZE-1);
	}
	int *ify = (int *)alloca((y1-y0+1) * sizeof(int));
	for (int y = y0; y <= y1; ++y) {
		real fy = abs((y - dImageY) * filter->invYWidth * FILTER_TABLE_SIZE);
		ify[y-y0] = MIN(Floor2Int(fy), FILTER_TABLE_SIZE-1);
	}

	for (int x = x0; x <= x1; ++x) {
		for (int y = y0; y <= y1; ++y)	{

			//////////////////////////////////////////////////////////////////////////
			// Evaluate filter value at $(x,y)$ pixel
			real filterWt = filterTable(ifx[x-x0], ify[y-y0]);

			//////////////////////////////////////////////////////////////////////////
			// Update pixel values with filtered sample contribution
			Pixel &pixel = pixels(x - xPixelStart, y - yPixelStart);

			pixel.L += (filterWt * L);
			pixel.alpha += alpha * filterWt;
			pixel.weightSum += filterWt;

			//nIter++;
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// Possibly write out in-progress image
	if (--sampleCount == 0) {
		writeImage(gamma);
		sampleCount = writeFrequency;
	}
}

void ImageFilm::getSampleExtent(int &xstart, int &xend, int &ystart, int &yend) const {
	xstart = Floor2Int(xPixelStart + .5f - filter->xWidth);
	xend   = Ceil2Int(xPixelStart - .5f + xPixelCount + filter->xWidth);
	ystart = Floor2Int(yPixelStart + .5f - filter->yWidth);
	yend   = Ceil2Int(yPixelStart - .5f + yPixelCount + filter->yWidth);
}

void ImageFilm::writeImage(real _gamma) {

	if(_gamma != 0.0f)
		gamma = _gamma;

	basis::Image image(xPixelCount, yPixelCount, basis::Image::IMAGE_RGB);

	RGBA color;
	real alpha;
	for(int x=0;x<xPixelCount;++x) {
		for(int y=0;y<yPixelCount;++y) {
			
			color[0] = pixels(x,y).L[0];	
			color[1] = pixels(x,y).L[1];	
			color[2] = pixels(x,y).L[2];	
			alpha = pixels(x,y).alpha;

			color /= pixels(x,y).weightSum;
			alpha /= pixels(x,y).weightSum;

			alpha = basis::clamp(alpha, static_cast<real>(0.0f), static_cast<real>(1.0f));
			if (premultiplyAlpha) 
				color *= alpha;

			color.clamp(0.0f, 1.0f);		

			image.setPixel(x,y,color);
		}
	}	

	if(gamma != 1.0f)
		image.gammaCorrection(gamma);
	image.saveAsTGA(filename.c_str());	
}

#ifdef _USE_XML_

void ImageFilm::importXmlNode(xmlNodePtr node) {
	Film::importXmlNode(node);

	xmlNodePtr filterNode = basis::XmlLibrary::xmlFindNode(node, "Filter");
	filter = Filter::createFromXmlNode(filterNode);

	xmlNodePtr cropNode = basis::XmlLibrary::xmlFindNode(node, "CropWindow");
	if(cropNode != NULL) {
		cropWindow[0] = basis::XmlLibrary::xmlGetFloat(cropNode, "Left");
		cropWindow[1] = basis::XmlLibrary::xmlGetFloat(cropNode, "Right");
		cropWindow[2] = basis::XmlLibrary::xmlGetFloat(cropNode, "Down");
		cropWindow[3] = basis::XmlLibrary::xmlGetFloat(cropNode, "Up");
	}
	else 
		throwError("No Crop Range _ illumination::ImageFilm::importXmlNode");

	filename = basis::XmlLibrary::xmlGetString(node, "Filename");
	premultiplyAlpha = basis::XmlLibrary::xmlGetBool(node, "PremultiplyAlpha");
	writeFrequency = sampleCount = basis::XmlLibrary::xmlGetInt(node, "WriteFrequency");

	init();
}

void ImageFilm::print_debug_info() const {
	const real zero(0.0);

	Spectrum _sum = Spectrum(zero);
	real _sum_alpha(zero), _sum_weight(zero);
	
	for(Array2<Pixel>::index_type i=0;i<pixels.dimX();++i) {
		for(Array2<Pixel>::index_type j=0;j<pixels.dimY();++j) {
			_sum += pixels(i,j).L;
			_sum_alpha += pixels(i,j).alpha;
			_sum_weight += pixels(i,j).weightSum;
		}
	}

	std::cout << "Spectrum Average _ " << _sum / (pixels.size()) << endl;
	std::cout << "Alpha Average _ " << _sum_alpha / (pixels.size()) << endl;
	std::cout << "Weight Average _ " << _sum_weight / (pixels.size()) << endl;
}

#endif