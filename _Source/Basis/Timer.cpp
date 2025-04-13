#include <stdafx.h>
#include <ctime>

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <mmsystem.h>
#include <Basis/Timer/Timer.h>

using namespace std;
using namespace basis;

void Timer::reset()	{
	startTime = timeGetTime() * 0.001f;
}

const double Timer::getTime()	const {
	if(isPaused)
		return pauseTime-startTime;
	else
		return (timeGetTime() * 0.001f) - startTime;
}

void Timer::printTime(std::string msg, std::ostream &os) {
	os << msg.c_str() << " : " << getTime() << "s" << std::endl;
}

void Timer::logTime(std::string msg, Log& _log) {
	_log << msg.c_str() << " : " << getTime() << "s" << std::endl;
}

void Timer::pause()	{
	if(isPaused)
		return;

	isPaused = true;
	pauseTime = timeGetTime() * 0.001f;
}

void Timer::unpause()	{
	if(!isPaused)
		return;

	isPaused = false;
	startTime += (timeGetTime() * 0.001f - pauseTime);
}

const std::string Timer::getCurrentTime() {
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

	return os.str();
}