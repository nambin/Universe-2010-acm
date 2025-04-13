#include <stdafx.h>
#include <Basis/TBBProgramming/TBBProgramming.h>

#ifdef USE_TBB

using namespace basis;
using namespace tbb;

TBBProgramming::TBBProgramming() {
	init = new task_scheduler_init(task_scheduler_init::deferred);
	nThreads = 0;
}

TBBProgramming::~TBBProgramming() {
	delete init;
	init = NULL;
}

void TBBProgramming::initialze(unsigned int _nThreads) {
	if(nThreads != 0)
		return;

	nThreads = _nThreads;
	init->initialize(nThreads);
}

void TBBProgramming::terminate() {
	if(nThreads == 0)
		return;

	init->terminate();
	nThreads = 0;
}

#endif		// USE_TBB