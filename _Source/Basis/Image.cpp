#include <stdafx.h>
#include <string>
#include <iostream>
#include <fstream>

#define _MFC_VER

extern "C" {
#include <libjpeg/jpeglib.h>
}

#include <Basis/Basis/BasicLibrary.h>
#include <Basis/StringLibrary/StringLibrary.h>
#include <Basis/DataPath/DataPath.h>
#include <Basis/Image/Image.h>

#ifdef _DEBUG
#define __ASSERT(c) assert(c)
#else
#define __ASSERT(c) {}
#endif

using namespace basis;

Image::Image() {
	pixels = NULL;
	width = height = 0;
	colorFormat = IMAGE_RGBA;
}

Image::Image(unsigned int w, unsigned int h, IMAGE_COLOR_FORMAT format) {
	if(w == 0 || h == 0) {
		pixels = NULL;
		width = height = 0;
		colorFormat = format;
		return;
	}

	width = w; height = h; colorFormat = format;
	pixels = new float[ colorFormat * width * height ];
	for ( unsigned int i = 0; i < colorFormat * width * height; ++i )
		pixels[i] = 0.0f;	
}

Image::~Image() {
	if(pixels)
		delete[] pixels;
}

void Image::resize(unsigned int w, unsigned int h, IMAGE_COLOR_FORMAT format) {
	if(pixels)
		delete[] pixels;

	if(w == 0 || h == 0) {
		pixels = NULL;
		width = height = 0;
		colorFormat = format;
		return;
	}

	width = w; height = h; colorFormat = format;
	pixels = new float[ colorFormat * width * height ];
	for ( unsigned int i = 0; i < colorFormat * width * height; ++i )
		pixels[i] = 0.0f;	
}

void Image::setPixel(unsigned int x, unsigned int y, const RGBA &color) {
	__ASSERT(x >= 0 && x < width && y >= 0 && y < height);

	for(unsigned int i=0;i<(unsigned int)colorFormat;++i)
		pixels[idx(x,y)*colorFormat + i] = color[i];
}

void Image::addPixel(unsigned int x, unsigned int y, const RGBA &color) {
	__ASSERT(x >= 0 && x < width && y >= 0 && y < height);

	for(unsigned int i=0;i<(unsigned int)colorFormat;++i)
		pixels[idx(x,y)*colorFormat + i] += color[i];
}

void Image::gammaCorrection(float gamma) {
	
	RGBA tmp;
	float power = 1.0f / gamma;

	for(unsigned int y=0;y<height;++y) 
		for(unsigned int x=0;x<width;++x) 
			for(unsigned int i=0;i<3;++i) 
				pixels[idx(x,y)*colorFormat + i] = pow(pixels[idx(x,y)*colorFormat + i], power);
}

bool Image::saveGLFrameBufferAsTGA( const char* const filename, unsigned int width, unsigned int height, unsigned int format) {

	unsigned int colorFormat;
	GLenum loadFormat;

	if ( format == GL_RGB || format == IMAGE_RGB ) {
		colorFormat = IMAGE_RGB;
		loadFormat = GL_RGB;
	}
	else if ( format == GL_RGBA || format == IMAGE_RGBA ) {
		colorFormat = IMAGE_RGBA;
		loadFormat = GL_RGBA;
	}
	else {
		basis::throwError("Unsupported Color Format\n");
		return 0;
	}	

	std::ofstream fout( filename, std::ios::binary );

	unsigned char TGAheader[12] = {0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0};	// Uncompressed TGA header
	unsigned char infoHeader[6];
	unsigned char *data = new unsigned char[colorFormat*width*height];

	glReadPixels( 0, 0, width, height, loadFormat, GL_UNSIGNED_BYTE, data );

	// Data needs to be in BGR format
	// Swap b and r
	for ( unsigned int i = 0; i < colorFormat*width*height; i += colorFormat ) {	
		// Repeated XOR to swap bytes 0 and 2
		data[i] ^= data[i+2] ^= data[i] ^= data[i+2];
	}

	// Save header
	for ( int i = 0; i < 12; i++ )
		fout << TGAheader[i];

	// Save info header
	infoHeader[0] = ( width  & 0x00FF );
	infoHeader[1] = ( width  & 0xFF00 ) >> 8;
	infoHeader[2] = ( height & 0x00FF );
	infoHeader[3] = ( height & 0xFF00 ) >> 8;
	infoHeader[4] = 8*colorFormat;
	infoHeader[5] = 0;

	for ( int i = 0; i < 6; i++ )
		fout << infoHeader[i];

	// write pixels
	fout.write(reinterpret_cast<const char*>(data), sizeof(unsigned char) * colorFormat * width * height);	
	fout.close();	

	delete[] data;

	return 1;
}

void Image::saveAsTGA(const char* const filename) const {

	std::ofstream fout( filename, std::ios::binary );

	unsigned char TGAheader[12] = {0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0};	//Uncompressed TGA header
	unsigned char infoHeader[6];

	unsigned char *data = new unsigned char[(int)colorFormat*width*height];

	for ( unsigned int j = 0; j < height; j++ ) {
		for ( unsigned int i = 0; i < width; i++ ) {		
			for( unsigned int c = 0;  c < (unsigned int)colorFormat;  c++ ) {
				float pixel = pixels[(unsigned int)colorFormat*idx(i,j)+c];

				// quantize
				pixel = pixel * 256.0f;
				pixel = pixel > 0.0f ? (pixel < 255.0f ? pixel : 255.0f) : 0.0f;

				data[(unsigned int)colorFormat*idx(i,j)+c] = (unsigned char)pixel;
			}
		}
	}

	// Data needs to be in BGR format
	// Swap b and r
	for ( unsigned int i = 0; i < (unsigned int)colorFormat*width*height; i += (unsigned int)colorFormat ) {	
		// Repeated XOR to swap bytes 0 and 2
		data[i] ^= data[i+2] ^= data[i] ^= data[i+2];
	}

	// Save header
	for ( int i = 0; i < 12; i++ )
		fout << TGAheader[i];

	// Save info header
	infoHeader[0] = ( width  & 0x00FF );
	infoHeader[1] = ( width  & 0xFF00 ) >> 8;
	infoHeader[2] = ( height & 0x00FF );
	infoHeader[3] = ( height & 0xFF00 ) >> 8;
	infoHeader[4] = 8*(int)colorFormat;
	infoHeader[5] = 8;
	
	for ( int i = 0; i < 6; i++ )
		fout << infoHeader[i];

	// write pixels
	fout.write(reinterpret_cast<const char*>(data), sizeof(unsigned char) * colorFormat * width * height);
	fout.close();

	delete[] data;
}

void Image::saveAsGLFrameBuffer() const {
	glDrawPixels( width, height, getGLColorFormat(), GL_FLOAT, pixels );
}

const int Image::load(const char* const filename) {

	std::string ext = basis::DataPath::getFilenameExtension(filename);

	if(ext.compare("jpg")==0 || ext.compare("jpeg")==0 
	|| ext.compare("JPG")==0 || ext.compare("JPEG")==0)		return loadFromJPEG(filename);
	else if(ext.compare("tga")==0 || ext.compare("TGA")==0)	return loadFromTGA(filename);

	return -1;
}

const int Image::loadFromTGA(const char* const filename) {

	std::ifstream fin(filename, std::ios::in | std::ios::binary );

	if ( !fin ) {
		std::cerr << "Error: Can't find file <" << filename << ">" << std::endl;
		return 0;
	}

	unsigned char header[12];
	unsigned char info[6];

	fin.read( (char*)header, 12 );
	fin.read( (char*)info, 6 );

	width = (unsigned int)info[0];
	width += ((unsigned int)info[1]) << 8;

	height = (unsigned int)info[2];
	height += ((unsigned int)info[3]) << 8;

	if ( (unsigned int)info[4] == 24 ) colorFormat = IMAGE_RGB;
	else if ( (unsigned int)info[4] == 32 ) colorFormat = IMAGE_RGBA;
	else {
		std::cerr << "Error: Unknown image colorFormat" << std::endl;
		return 0;
	}

	unsigned int format = (unsigned int)colorFormat;

	// Create new raster with correct size 
	if ( pixels ) delete [] pixels;
	pixels = new float[format*width*height];
	char *data = new char[format*width*height];

	fin.read( data, format*width*height );

	// Data needs to be in BGR format
	// Swap b and r
	for ( unsigned int i = 0; i < format*width*height; i += format ) {	
		// Repeated XOR to swap bytes 0 and 2
		data[i] ^= data[i+2] ^= data[i] ^= data[i+2];
	}

	for ( unsigned int i = 0; i < width; i++ ) {
		for ( unsigned int j = 0; j < height; j++ ) {
			for( unsigned int c = 0;  c < format; c++ ) {
				pixels[format*idx(i,j)+c] = (float)((unsigned char)data[format*idx(i,j)+c])/256.f;
			}
		}
	}

	fin.close();
	delete [] data;

	return 1;
}

const int Image::loadFromJPEG( const char *filename ) {
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;

	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_decompress(&cinfo);

	FILE* fp;
	fopen_s(&fp, filename, "rb");

	if ( fp == NULL ) {
		std::cerr << "Error: Can't find file <" << filename << ">" << std::endl;
		return 0;
	}

	jpeg_stdio_src(&cinfo, fp);
	jpeg_read_header(&cinfo, true);

	width = cinfo.image_width;
	height = cinfo.image_height;
	int d = cinfo.jpeg_color_space;

	jpeg_start_decompress(&cinfo);
	unsigned char *data = new unsigned char[width * height * d];
	while (cinfo.output_scanline < cinfo.output_height)	{
		jpeg_read_scanlines(&cinfo, &data, 1);
		data += d * cinfo.output_width;
	}
	data -= d * cinfo.output_width * cinfo.output_height;

	if ( pixels ) delete [] pixels;
	pixels = new float[width * height * d];
	for ( unsigned int i = 0; i < width; i++ ) {
		for ( unsigned int j = 0; j < height; j++ ) {
			for( int c = 0;  c < d; c++ ) {
				pixels[d*idx(i,j)+c] = (float)((unsigned char)data[d*idx(i,height-j-1)+c])/256.f;
			}
		}
	}

	if ( d == 3 ) colorFormat = IMAGE_RGB;
	else if ( d == 4 ) colorFormat = IMAGE_RGBA;

	fclose(fp);
	jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);
	delete [] data;

	return 1;
}

const int Image::loadFromGLFrameBuffer(int x, int y, int w, int h, unsigned int format) {

	if( w <= 0 && h <= 0 ) {
		basis::throwError("Wrong Image Size\n");
		return 0;
	}
	
	width = w;
	height = h;

	GLenum loadFormat;

	if ( format == GL_RGB || format == IMAGE_RGB ) {
		colorFormat = IMAGE_RGB;
		loadFormat = GL_RGB;
	}
	else if ( format == GL_RGBA || format == IMAGE_RGBA ) {
		colorFormat = IMAGE_RGBA;
		loadFormat = GL_RGBA;
	}
	else {
		basis::throwError("Unsupported Color Format\n");
		return 0;
	}		

	if ( pixels ) delete [] pixels;
	pixels = new float[ (unsigned int)colorFormat * width * height ];

	glReadPixels( x, y, width, height, loadFormat, GL_FLOAT, pixels );

	return 1;
}

#undef __ASSERT