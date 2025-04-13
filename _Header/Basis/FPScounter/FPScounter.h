#ifndef _BASIS_FPS_COUNTER_H_
#define _BASIS_FPS_COUNTER_H_

#include "../Basis.h"

namespace basis {

class BASISDLL FPScounter	{
public:
	FPScounter() : fps(0.0f), lastTime(0.0f), frames(0), time(0.0f) {}

	void update();
	const double getFps() const { return fps; }

protected:
	double fps;

	double lastTime;	
	double time;

	int frames;
};

};	// namespace basis

#endif	