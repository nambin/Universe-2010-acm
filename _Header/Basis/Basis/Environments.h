#pragma once

#include <string>
#include "../Basis.h"
#include "Log.h"

//#define USE_TBB
//#define USE_OPENMP

#ifdef _OPENMP
#undef USE_OPENMP
#define USE_OPENMP
#endif

#ifdef USE_OPENMP
#include <omp.h>

#define NUM_OPENMP_THREADS_MAX 8
#endif

static const std::string output_folder = "../Output/";

BASISDLL extern basis::Log log_info;