#ifndef _BASIS_IMAGE_H_
#define _BASIS_IMAGE_H_

#include <cstdlib>
#ifdef _WIN32
#include "gl/glut.h"
#else
#include <GLUT/glut.h>
#endif

#include "../Color/RGBA.h"

namespace basis {

//////////////////////////////////////////////////////////////////////////
// Image Class Definition
//////////////////////////////////////////////////////////////////////////
class BASISDLL Image {
public:
	enum IMAGE_COLOR_FORMAT {
		IMAGE_RGB = 3, 
		IMAGE_RGBA = 4
	};

	////////////////////////////////////////////////////////////////////////////
	// Constructors / Destructors
	Image();
	Image( unsigned int w, unsigned int h, IMAGE_COLOR_FORMAT format = IMAGE_RGBA );
	~Image();

	////////////////////////////////////////////////////////////////////////////
	// Setter
	void resize( unsigned int w, unsigned int h, IMAGE_COLOR_FORMAT format = IMAGE_RGBA );

	void setPixel( unsigned int x, unsigned int y, const RGBA & );
	void addPixel( unsigned int x, unsigned int y, const RGBA & );

	////////////////////////////////////////////////////////////////////////////
	// Getter
	const unsigned int getWidth() const { return width; }
	const unsigned int getHeight() const { return height; }
	const IMAGE_COLOR_FORMAT getColorFormat() const { return colorFormat; }
	const GLenum getGLColorFormat() const {
		if ( colorFormat == IMAGE_RGBA ) return GL_RGBA;
		else return GL_RGB;
	}

	////////////////////////////////////////////////////////////////////////////
	// Manipulation
	void gammaCorrection(float gamma);

	////////////////////////////////////////////////////////////////////////////
	// Get Pointer
	const float* const getPixels() const { return pixels; }
	const float* const getptr() const { return pixels; }
	
	////////////////////////////////////////////////////////////////////////////
	// Save	
	void saveAsTGA( const char* const filename ) const;
	void saveAsGLFrameBuffer() const;
	
	////////////////////////////////////////////////////////////////////////////
	// Load
	const int load(const char* const filename);
	const int loadFromTGA( const char* const filename );
	const int loadFromJPEG( const char* const filename );
	const int loadFromGLFrameBuffer( int x, int y, int w, int h, unsigned int format );

	////////////////////////////////////////////////////////////////////////////
	// Static Function
	static bool saveGLFrameBufferAsTGA( const char* const filename, unsigned int width, unsigned int height, unsigned int format);

private : 
	const unsigned int idx(unsigned int x, unsigned int y) const { return y * width + x; }

private:
	float					*pixels;		
	unsigned int			width;			
	unsigned int			height;	

	IMAGE_COLOR_FORMAT		colorFormat;
};

};		// namespace basis

#endif 
