#pragma once

#include "../Basis/Environments.h"

#ifdef USE_TBB

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#include <tbb/blocked_range2d.h>
#include <tbb/blocked_range3d.h>
#include <tbb/tick_count.h>
#include "../Basis.h"

namespace basis {

class BASISDLL TBBProgramming {
public :
	TBBProgramming();
	~TBBProgramming();

	void initialze(unsigned int nThreads);
	void terminate();

public : 
	tbb::task_scheduler_init*	init;
	unsigned int				nThreads;
};

};		// basis

#endif		// USE_TBB