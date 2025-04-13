#include <stdafx.h>
#include <fstream>
#include <sstream>
#include <ctime>
#include <Basis/Basis/Log.h>
#include <Basis/Basis/Environments.h>

using namespace std;
using namespace basis; 

basis::Log log_info("log", output_folder.c_str(), true);

Log::Log(const string& file, const string& fold, const bool _is_print, unsigned int freq) : filename(file), folder(fold), is_print(_is_print), num_acc_logs(0), freq_flush(freq), ifs("") {
	make_full_filename();
	ifs.close();
	ifs.clear();
}

void Log::set_folder(const string& fold) {
	if(folder.compare(fold) == 0) 
		return;

	folder = fold;
	make_full_filename();

	if(ifs.is_open() == true)
		make_fileIO();
}

void Log::set_filename(const string& file) {
	if(filename.compare(file) == 0)
		return;

	filename = file;
	make_full_filename();

	if(ifs.is_open() == true)
		make_fileIO();
}

void Log::make_full_filename() {
	time_t timer; 
	struct tm *t = new struct tm();

	timer = time(NULL);
	localtime_s(t, &timer); 	

	std::stringstream os;
	os << (t->tm_year+1900);	os << "-";
	os << (t->tm_mon+1);		os << "-";
	os << t->tm_mday;			os << " ";

	os << t->tm_hour;			os << "-";
	os << t->tm_min;			os << "";

	filename_full = folder + filename + " _ " + os.str() + ".txt";

	delete t;
}

void Log::make_fileIO() {
	ifs.close(); 
	ifs.clear();
	ifs.open(filename_full.c_str(), fstream::in | ios::app);
}

void Log::add_message(StandardEndLine manip)    {			
	if(ifs.is_open() == false)
		make_fileIO();

	ifs << os.str(); os.str(""); os.clear();
	ifs << "\n"; num_acc_logs++;

	if(is_print)
		manip(std::cout);	

	if(num_acc_logs % freq_flush == 0)
		ifs.close();
}