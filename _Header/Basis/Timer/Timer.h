#ifndef _BASIS_TIMER_H_
#define _BASIS_TIMER_H_

#include <iostream>
#include "../Basis.h"
#include "../Basis/Log.h"

namespace basis {

class BASISDLL Timer	{
public:
	Timer() : isPaused(false) { reset(); }

	void reset();
	const double getTime() const;

	void printTime(std::string msg, std::ostream &os = std::cout);
	void logTime(std::string msg, Log& _log);

	void pause();
	void unpause();

	static const std::string getCurrentTime();

protected:
	bool isPaused;

	double startTime;	
	double pauseTime;
};

}

#endif	