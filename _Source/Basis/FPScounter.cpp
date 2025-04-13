#include <stdafx.h>

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <mmsystem.h>
#include <Basis/FPScounter/FPScounter.h>

using namespace std;
using namespace basis;

void FPScounter::update()	{

	// keep track of time passed and frame count
	time = timeGetTime() * 0.001f;
	++frames;

	//If a second has passed
	if(time-lastTime > 1.0f) {
		fps = frames / (time-lastTime);	//update the number of frames per second
		lastTime = time;				//set time for the start of the next count
		frames = 0;						//reset fps for this second
	}
}