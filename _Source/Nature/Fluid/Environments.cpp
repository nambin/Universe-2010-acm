#include <stdafx.h>
#include <Nature/Nature.h>
#include <Nature/Fluid/Environment/Environments.h>

using namespace std;
using namespace basis;

namespace fluidspace {

const std::string get_concat_filename(const std::string& name, const std::size_t frame) {
	char filename[512];
	sprintf_s(filename, std::string(" _ %04d").c_str(), frame);

	return name + std::string(filename);
}

};