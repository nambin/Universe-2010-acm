#pragma once

#include <iostream>
#include <string>
#include <algorithm>							
#include "TBBProgramming/TBBProgramming.h"

#ifdef USE_TBB

using namespace tbb;
using namespace std;

namespace basis {

BASISDLL void TestSubStringFinder();
BASISDLL void TestMatrixMultiply();

inline void TestTBB() {
//	TestSubStringFinder();
//	TestMatrixMultiply();	
}

};		// basis

#endif		// USE_TBB