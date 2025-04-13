#ifndef _SHADERSPACE_CGSHADER_VOLUME_RENDERING_H_
#define _SHADERSPACE_CGSHADER_VOLUME_RENDERING_H_

#include "../../IlluminationSpace.h"
#include <cstdlib>
#include <iostream>
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glut.h>
#include "../../../Basis/Color/RGBA.h"
#include "../../../Principle/Arithmetic/Vector/Vector3.h"
#include "../../../Principle/Arithmetic/GridLibrary/Grid3.h"
#include "../../../Principle/Geometry/DistanceFieldLibrary/DistanceFields.h"
#include "../../../Principle/Geometry/MacGrid/MacGridLibrary.h"
#include "../../../Principle/VirtualSpace/DrawLibrary/glDrawLibrary.h"
#include "../../../Principle/VirtualSpace/Scene/Scene.h"
#include "cgShaderLibrary.h"

namespace shaderspace {

using basis::RGBA;
using arithmetic::Vector3;
using geometryspace::DistanceField3D;

class ILLUMINATIONDLL cgShader_VolumeRendering {
public : 
	typedef arithmetic::Grid3<float,float>::index_type index_type;

public : 
	cgShader_VolumeRendering();

	template<typename T>
	void setVolumeData(const arithmetic::Grid3<T,T>&, const RGBA& c = RGBA::white());

	template<typename T>
	void setVolumeData(const Vector3<index_type>& width, const Vector3<T>& L, const DistanceField3D<T>&, const RGBA& c = RGBA::white());

	void glDraw() const;
	void resizeScreen(unsigned int width, unsigned int height);

	// Auxiliary Function
	void setVolumeData_Default();

private : 
	void init();
	void glewInit();
	void cgInit();
	void FBOInit(unsigned int width, unsigned int height);	

	void createVolumeTexture(GLsizei width, GLsizei height, GLsizei depth, GLubyte* data);

	void enableRenderbuffers() const;
	void disableRenderBuffers() const;

	void renderBackFace() const;
	void rayCasting() const;

	void drawVertex(float x, float y, float z) const;
	void drawQuads(float x, float y, float z) const;
	
public : 
	bool					dBoundingBox;

	CGcontext				context; 
	CGprofile				vertexProfile, fragmentProfile; 
	CGprogram				vertexProgram, fragmentProgram;
	mutable CGparameter		param1,param2;

	GLuint					volume_texture;		// the volume texture	
	float					Lx, Ly, Lz;
	float					stepsize;

	GLuint					frameBuffer; 
	GLuint					backFaceBuffer;		// the FBO buffers	
};

template<typename T>
void cgShader_VolumeRendering::setVolumeData(const arithmetic::Grid3<T,T>& grid, const RGBA& _c) { 	

	index_type width = grid.nX;
	index_type height = grid.nY;
	index_type depth = grid.nZ;

	Lx = grid.dx * (width-1); 
	Ly = grid.dy * (height-1); 
	Lz = grid.dz * (depth-1);

	index_type size = width * height * depth * 4;
	GLubyte *data = new GLubyte[size];

	RGBA c = _c;
	c[0] = MIN(c[0]*255.0f, float(255)); c[1] = MIN(c[1]*255.0f, float(255));
	c[2] = MIN(c[2]*255.0f, float(255)); c[3] = MIN(c[3]*255.0f, float(255));

	index_type idx;
	for(index_type x = 0; x < width; x++) {
		for(index_type y = 0; y < height; y++) {
			for(index_type z = 0; z < depth; z++) {
				idx = ((x * height + y) * depth) + z;
				data[idx*4+0] = MIN(c[0], float(255));
				data[idx*4+1] = MIN(c[1], float(255));
				data[idx*4+2] = MIN(c[2], float(255));
				data[idx*4+3] = MIN(float(grid(x,y,z) * c[3]), float(255));
			}
		}
	}

	createVolumeTexture(static_cast<GLsizei>(width), static_cast<GLsizei>(height), static_cast<GLsizei>(depth), data);

	delete [] data;
}

template<typename T>
void cgShader_VolumeRendering::setVolumeData(const Vector3<index_type>& _width, const Vector3<T>& L, const DistanceField3D<T>& field, const RGBA& _c) {

	index_type width = _width[0];
	index_type height = _width[1];
	index_type depth = _width[2];

	Lx = L[0]; Ly = L[1]; Lz = L[2];

	float dx = Lx / (width-1);
	float dy = Ly / (height-1);
	float dz = Lz / (depth-1);
	float hInv = MIN3(1.0f/dx, 1.0f/dy, 1.0f/dz);

	index_type size = width * height * depth * 4;
	GLubyte *data = new GLubyte[size];

	RGBA c = _c;
	c[0] = MIN(c[0]*255.0f, float(255)); c[1] = MIN(c[1]*255.0f, float(255));
	c[2] = MIN(c[2]*255.0f, float(255)); c[3] = MIN(c[3]*255.0f, float(255));

	float alpha, dist;
	index_type idx;
	Vector3<float> pos;
	for(index_type x = 0; x < width; x++) {
		for(index_type y = 0; y < height; y++) {
			for(index_type z = 0; z < depth; z++) {
				pos = Vector3<T>(dx * x, dy * y, dz * z);
				dist = field.getPhi(pos);
				alpha = MacGridLibrary<T>::heaviside(dist * hInv, T(0), T(1));	

				idx = ((x * height + y) * depth) + z;
				data[idx*4+0] = MIN(c[0], float(255));
				data[idx*4+1] = MIN(c[1], float(255));
				data[idx*4+2] = MIN(c[2], float(255));
				data[idx*4+3] = MIN(alpha * c[3], float(255));
			}
		}
	}

	createVolumeTexture(static_cast<GLsizei>(width), static_cast<GLsizei>(height), static_cast<GLsizei>(depth), data);

	delete [] data;
}

};			// shaderspace

#endif