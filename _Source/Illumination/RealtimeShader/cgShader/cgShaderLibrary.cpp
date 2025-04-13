#include <stdafx.h>
#include <Illumination/RealtimeShader/cgShader/cgShaderLibrary.h>

using namespace std;
using namespace shaderspace;

bool cgShaderLibrary::set_cgProfile(CGprofile& p, const CGprofile profile) {
	if (cgGLIsProfileSupported(profile)) {
		p = profile;
		return true;		
	}
	return false;
}

bool cgShaderLibrary::load_cgProgram(CGprogram& program, CGcontext context, CGprofile profile, string shader_path, string program_name) {
	assert(cgIsContext(context));

	program = cgCreateProgramFromFile(	context, 
										CG_SOURCE,
										shader_path.c_str(),
										profile,
										program_name.c_str(), 
										NULL);	

	if (!cgIsProgramCompiled(program))
		cgCompileProgram(program);

	cgGLEnableProfile(profile);
	cgGLLoadProgram(program);
	cgGLDisableProfile(profile);

	return true;
}

void cgShaderLibrary::set_TextureParameter(CGparameter& param, string str, GLuint tex, const CGprogram &program) {
	param = cgGetNamedParameter(program, str.c_str()); 
	cgGLSetTextureParameter(param, tex); 
	cgGLEnableTextureParameter(param);
}