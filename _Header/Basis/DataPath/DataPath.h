#ifndef _BASIS_DATAPATH_H_
#define _BASIS_DATAPATH_H_

#ifdef _WIN32
#pragma warning(disable:4786)   // symbol size limitation ... STL
#endif

#include <vector>
#include <string>
#include "../Basis.h"

namespace basis {

class BASISDLL DataPath	{
public:
	DataPath() {};
	~DataPath() {};

	void putPath(const std::string path);
	const std::string getPath(const std::string filename);
	const std::string getFile(const std::string filename);

	static const std::string getDirectory(const std::string filename);
	static const std::string getFilenameExtension(const std::string filename);

public : 
	std::string              file_path;
	std::string              path_name;
	std::vector<std::string> path;
};

};		// namespace basis

#endif
