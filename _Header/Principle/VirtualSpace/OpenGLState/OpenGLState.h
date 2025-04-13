#ifndef _VIRTUALSPACE_OPENGL_STATE_H_
#define _VIRTUALSPACE_OPENGL_STATE_H_

#include "OpenGLStateVariable.h"
#include "OpenGLBinaryStateVariable.h"
#include "../../Principle.h"
#include <vector>

namespace virtualspace {

class PRINCIPLEDLL OpenGLState {
public : 
	static const int numBinaryState = 8;
	enum IDX_BINARY_STATE_VARIABLE {
		IDX_BINARY_LIGHTING = 0, IDX_BINARY_DEPTHTEST, IDX_BINARY_CULLFACE, IDX_BINARY_BLEND, IDX_BINARY_POINTSMOOTH, 
		IDX_BINARY_LINESMOOTH, IDX_BINARY_POLYGONSMOOTH, IDX_BINARY_COLORMATERIAL
	};

	static const int numTwoInputState = 3;
	enum IDX_TWOINPUT_STATE_VARIABLE {
		IDX_TWOINPUT_COLORMATERIAL = 0, IDX_TWOINPUT_POLYGONMODE, IDX_TWOINPUT_BLENDFUNC
	};

public : 
	OpenGLState();
	OpenGLState(const OpenGLState &);
	void operator=(const OpenGLState &);

	void init();
	void glBind() const;

	void enable_Lighting()			{ binaryState[IDX_BINARY_LIGHTING]->enable = true; }
	void disable_Lighting()			{ binaryState[IDX_BINARY_LIGHTING]->enable = false; }

	void enable_DepthTest()			{ binaryState[IDX_BINARY_DEPTHTEST]->enable = true; }
	void disable_DepthTest()		{ binaryState[IDX_BINARY_DEPTHTEST]->enable = false; }

	void enable_CullFace()			{ binaryState[IDX_BINARY_CULLFACE]->enable = true; }
	void disable_CullFace()			{ binaryState[IDX_BINARY_CULLFACE]->enable = false; }

	void enable_Blend()				{ binaryState[IDX_BINARY_BLEND]->enable = true; }
	void disable_Blend()			{ binaryState[IDX_BINARY_BLEND]->enable = false; }

	void enable_PointSmooth()		{ binaryState[IDX_BINARY_POINTSMOOTH]->enable = true; }
	void disable_PointSmooth()		{ binaryState[IDX_BINARY_POINTSMOOTH]->enable = false; }

	void enable_LineSmooth()		{ binaryState[IDX_BINARY_LINESMOOTH]->enable = true; }
	void disable_LineSmooth()		{ binaryState[IDX_BINARY_LINESMOOTH]->enable = false; }

	void enable_PolygonSmooth()		{ binaryState[IDX_BINARY_POLYGONSMOOTH]->enable = true; }
	void disable_PolygonSmooth()	{ binaryState[IDX_BINARY_POLYGONSMOOTH]->enable = false; }

	void enable_ColorMaterial()		{ binaryState[IDX_BINARY_COLORMATERIAL]->enable = true; }
	void disable_ColorMaterial()	{ binaryState[IDX_BINARY_COLORMATERIAL]->enable = false; }

	void set_ColorMaterial(GLenum face, GLenum mode)	{ twoInputState[IDX_TWOINPUT_COLORMATERIAL]->set(face, mode);	}
	void set_PolygonMode(GLenum face, GLenum mode)		{ twoInputState[IDX_TWOINPUT_POLYGONMODE]->set(face, mode);		}
	void set_BlendFunc(GLenum s, GLenum d)				{ twoInputState[IDX_TWOINPUT_BLENDFUNC]->set(s,d);				}

private : 
	std::vector<OpenGLBinaryStateVariable*>		binaryState;
	std::vector<OpenGLTwoInputStateVariable*>	twoInputState;
};

};		// namespace virtualspace

#endif