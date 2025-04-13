#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include "../Basis.h"

namespace basis {

class BASISDLL Log {
public : 
	typedef std::basic_ostream<char, std::char_traits<char>> CoutType;		// this is the function signature of std::endl    
	typedef CoutType& (*StandardEndLine)(CoutType&);						// define an operator<< to take in std::endl   

	Log(const std::string& filename, const std::string& folder = "./", const bool is_print = false, unsigned int freq_flush = 250);

	void set_folder(const std::string&);
	void set_filename(const std::string&);
	void set_flush_frequency(unsigned int freq) { freq_flush = freq; }
	
	void enable_print() { is_print = true; }
	void disable_print() { is_print = false; }

	template<typename MSG>
	void add_message(const MSG&);

	void add_message(StandardEndLine manip);

private : 
	void make_fileIO();	
	void make_full_filename();

private : 
	std::string				folder;
	std::string				filename;

	bool					is_print;
	unsigned int			freq_flush;

	std::fstream			ifs;
	std::string				filename_full;

	std::stringstream		os;
	mutable unsigned int	num_acc_logs;	
};

template<typename MSG>
Log& operator<<(Log &_log, const MSG& msg) {
	_log.add_message(msg);
	return _log;
}

inline Log& operator<<(Log &_log, const Log::StandardEndLine msg) {
	_log.add_message(msg);
	return _log;
}

template<typename MSG>
void Log::add_message(const MSG& msg) {
	if(ifs.is_open() == false)
		make_fileIO();	

	os << msg;	num_acc_logs++;

	if(is_print)
		std::cout << msg;

	if(num_acc_logs % freq_flush == 0)
		ifs.close();
}

};			// basis