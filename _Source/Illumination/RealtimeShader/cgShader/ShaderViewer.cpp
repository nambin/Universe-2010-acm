#include <stdafx.h>
#include <Illumination/RealtimeShader/ShaderViewer/ShaderViewer.h>
#include <Principle/Geometry/DistanceFieldLibrary/DistanceFields.h>
#include <Principle/VirtualSpace/Scene/Scene.h>

using namespace shaderspace;
using namespace geometryspace;
using namespace virtualspace;

ShaderViewer::ShaderViewer() {
	ioHandler = new ShaderViewer_IOHandler(this);
	//vr.setVolumeData_Default();
	//vr.setVolumeData(Vector3<size_t>(100,200,100), Vector3<float>(1,2,1), SphereDistanceField3D<float>(Vector3<float>(0.5f,0.3f,0.5f), 0.2f), RGBA(1,0,0,0.6f));

	SphereDistanceField3D<float> field(Vector3<float>(0.5f,0.7f,0.5f), float(0.2f));
	vr.setVolumeData(Vector3<cgShader_VolumeRendering::index_type>(128), Vector3<float>(1), field);
}

ShaderViewer::~ShaderViewer() {
}

void ShaderViewer::glBind() const {
	glDraw();
}

void ShaderViewer::glDraw() const {
	vr.glDraw();
}

void ShaderViewer::resizeScreen(unsigned int _width, unsigned int _height) {
	vr.resizeScreen(_width, _height);
}