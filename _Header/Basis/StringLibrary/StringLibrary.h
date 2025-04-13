#ifndef _BASIS_STRING_LIBRARY_H_
#define _BASIS_STRING_LIBRARY_H_

#include <string>
#include <atlconv.h>
#include "../Basis.h"

namespace basis {

class BASISDLL StringLibrary {
public : 
	static void replace(std::string& str, const std::string &src, const std::string &dest);
	static const size_t count(std::string& str, const std::string &src);

	static const std::wstring string2wstring(const std::string& src);
	static const std::string wstring2string(const std::wstring& src);

	static const std::string base64_encode(const std::string& src);
	static const std::string base64_decode(const std::string& src);

	static void importString_Binary(std::string&, std::ifstream& ifs);
	static void exportString_Binary(const std::string&, std::ofstream& ofs);
};

};		// namespace basis

#endif

