#ifndef _SHADERSPACE_CG_SHADER_LIBRARY_H_
#define _SHADERSPACE_CG_SHADER_LIBRARY_H_

#include "../../IlluminationSpace.h"
#include <string>
#include <iostream>
#include <cassert>
#include <Cg/cg.h>
#include <Cg/cgGL.h>

namespace shaderspace {

using std::string;

class ILLUMINATIONDLL cgShaderLibrary {
public : 
	static bool set_cgProfile(CGprofile&, const CGprofile);
	static bool load_cgProgram(CGprogram& program, CGcontext context, CGprofile profile, string shader_path, string program_name);
	static void set_TextureParameter(CGparameter& param, string str, GLuint tex, const CGprogram &program);
};

};		// shaderspace

#endif