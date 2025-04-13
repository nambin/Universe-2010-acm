#include <stdafx.h>
#include <Illumination/RealtimeShader/cgShader/cgShader_VolumeRendering.h>
#include <Principle/VirtualSpace/Scene/Scene.h>

using namespace std;
using namespace basis;
using namespace arithmetic;
using namespace shaderspace;
using namespace virtualspace;

cgShader_VolumeRendering::cgShader_VolumeRendering() { 
	backFaceBuffer = frameBuffer = volume_texture = -1;
	dBoundingBox = true;
	init(); 
}

void cgShader_VolumeRendering::init() {
	Lx = Ly = Lz = 1.0f;
	stepsize = 1.0f / 50.0f;

	glewInit();
	cgInit();		
}

void cgShader_VolumeRendering::glewInit() {
	GLenum err = ::glewInit();

	// initialize all the OpenGL extensions
	glewGetExtension("glMultiTexCoord2fvARB");  

	if(glewGetExtension("GL_EXT_framebuffer_object"))	cout << "GL_EXT_framebuffer_object support " << endl;
	if(glewGetExtension("GL_EXT_renderbuffer_object"))	cout << "GL_EXT_renderbuffer_object support " << endl;
	if(glewGetExtension("GL_ARB_vertex_buffer_object")) cout << "GL_ARB_vertex_buffer_object support" << endl;
	if(GL_ARB_multitexture)								cout << "GL_ARB_multitexture support " << endl;

	if (glewGetExtension("GL_ARB_fragment_shader")      != GL_TRUE ||
		glewGetExtension("GL_ARB_vertex_shader")        != GL_TRUE ||
		glewGetExtension("GL_ARB_shader_objects")       != GL_TRUE ||
		glewGetExtension("GL_ARB_shading_language_100") != GL_TRUE)
	{
		throwError("Driver does not support OpenGL Shading Language\n");
	}
}

void cgShader_VolumeRendering::cgInit() {
	context = cgCreateContext();

	if(cgShaderLibrary::set_cgProfile(vertexProfile, CG_PROFILE_VP40))			cout << "CG_PROFILE_VP40 supported." << endl; 
	else if(cgShaderLibrary::set_cgProfile(vertexProfile, CG_PROFILE_ARBVP1))	cout << "CG_PROFILE_ARBVP1 supported." << endl; 
	else 
		throwError("Neither arbvp1 or vp40 vertex profiles supported on this system\n");

	if(cgShaderLibrary::set_cgProfile(fragmentProfile, CG_PROFILE_FP40))		cout << "CG_PROFILE_FP40 supported." << endl;
	else if(cgShaderLibrary::set_cgProfile(fragmentProfile, CG_PROFILE_ARBFP1))	cout << "CG_PROFILE_ARBFP1 supported." << endl;
	else 
		throwError("Neither arbfp1 or fp40 fragment profiles supported on this system\n");

	// load the vertex and fragment raycasting programs
	cgShaderLibrary::load_cgProgram(vertexProgram, context, vertexProfile, "../Illumination/cgShader/volumeShader.cg", "vertex_main");
	cgShaderLibrary::load_cgProgram(fragmentProgram, context, fragmentProfile, "../Illumination/cgShader/volumeShader.cg", "fragment_main");
}

void cgShader_VolumeRendering::resizeScreen(unsigned int width, unsigned int height) {
	FBOInit(width, height);
}
	
void cgShader_VolumeRendering::FBOInit(unsigned int width, unsigned int height) {
	glGenFramebuffersEXT(1, &frameBuffer);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, frameBuffer);

	if(backFaceBuffer != -1)
		glDeleteTextures(1, &backFaceBuffer);

	glGenTextures(1, &backFaceBuffer);
	glBindTexture(GL_TEXTURE_2D, backFaceBuffer);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F_ARB, width, height, 0, GL_RGBA, GL_FLOAT, NULL);

//	cout << backFaceBuffer << endl;
}
