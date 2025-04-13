#include <stdafx.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <Basis/DataPath/DataPath.h>

using namespace std;
using namespace basis;

void DataPath::putPath(const std::string _path) {
	path.push_back(_path);
}

const std::string DataPath::getPath(const std::string filename) {

	FILE* fp;
	bool found = false;
	for(unsigned int i=0; i < path.size(); ++i)	{
		path_name = path[i] + "/" + filename;
		::fopen_s(&fp, path_name.c_str(), "r");

		if(fp != 0)	{     
			fclose(fp);
			found = true;
			break;
		}
	}

	if (found == false)	{
		path_name = filename;
		::fopen_s(&fp, path_name.c_str(), "r");
		if (fp != 0) {
			fclose(fp);
			found = true;
		}
	}

	if (found == false)	{
		for(unsigned int i=0; i < path.size(); ++i)	{
			unsigned int pos; 
			string fname = filename;
			while ((pos = (unsigned int)fname.find_first_of("\\/")) != string::npos)	{
				pos++; 

				fname = fname.substr(pos, fname.length()-pos);
				path_name = path[i] + "/" + fname;
				::fopen_s(&fp, path_name.c_str(), "r");

				if(fp != 0)	{     
					fclose(fp);
					found = true;
					break;
				}
			}

			if (found) break;
		}
	}

	if (found == false)
		return "";

	int loc = (int)path_name.rfind('\\');
	if (loc == -1)
		loc = (int)path_name.rfind('/');

	if (loc != -1)
		file_path = path_name.substr(0, loc);
	else
		file_path = ".";

	return file_path;
}

const std::string DataPath::getFile(const std::string filename) {

	FILE* fp;

	for(unsigned int i=0; i < path.size(); ++i)	{
		path_name = path[i] + "/" + filename;
		::fopen_s(&fp, path_name.c_str(), "r");

		if(fp != 0)	{     
			fclose(fp);
			return path_name;
		}
	}

	path_name = filename;
	::fopen_s(&fp, path_name.c_str(),"r"); 
	if (fp != 0) {
		fclose(fp);
		return path_name;
	}

	for(unsigned int i=0; i < path.size(); ++i)	{
		unsigned int pos; 
		std::string fname = filename;
		while ((pos = (unsigned int)fname.find_first_of("\\/")) != std::string::npos) {
			pos++; 

			fname = fname.substr(pos, fname.length()-pos);
			path_name = path[i] + "/" + fname;
			::fopen_s(&fp, path_name.c_str(), "r");

			if(fp != 0) {     
				fclose(fp);
				return path_name;
			}
		}
	}

	std::cerr << "File not found : " << filename << std::endl;
	return "";
}

const std::string DataPath::getDirectory(const std::string filename) {

	std::string ret(filename);

	size_t found = ret.rfind("/");
	if(found == string::npos)
		found = ret.rfind("\\");

	if(found != string::npos)
		ret.replace(found+1,ret.size()-found-1, "");
	else
		ret.replace(0,ret.size(), "");

	return ret;
}

const std::string DataPath::getFilenameExtension(const std::string filename) {

	std::string ret(filename);

	size_t found = ret.rfind(".");
	if(found == string::npos)
		return std::string("");

	ret.replace(0, found+1, "");	

	return ret;

}