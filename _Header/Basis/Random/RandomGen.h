#ifndef _BASIS_RANDOM_H_
#define _BASIS_RANDOM_H_

#include "../Basis/Environments.h"
#include "../Randomc/randomc.h"

namespace basis {

template<typename T>
class RandomGen {
public:
	RandomGen(int seed = 0) { setSeed(seed); }

#ifdef USE_OPENMP
	void setSeed(int seed) { int seeds[2] = { seed, 0 }; for(int i=0;i<NUM_OPENMP_THREADS_MAX;++i) { seeds[1] = i; mersenne[i].RandomInitByArray(seeds, 2); } }

	const int getInt() const			{ return mersenne[omp_get_thread_num()].IRandom(INT_MIN, INT_MAX);		}
	const float getFloat() const		{ return static_cast<float>(mersenne[omp_get_thread_num()].Random());	}
	const double getDouble() const		{ return mersenne[omp_get_thread_num()].Random();						}

	const T operator()() const			{ return static_cast<T>(mersenne[omp_get_thread_num()].Random());		}

private:
	mutable CRandomMersenne mersenne[NUM_OPENMP_THREADS_MAX];

#else
	void setSeed(int seed) { mersenne.RandomInit(seed); }

	const int getInt() const			{ return mersenne.IRandom(INT_MIN, INT_MAX);	}
	const float getFloat() const		{ return static_cast<float>(mersenne.Random());	}
	const double getDouble() const		{ return mersenne.Random();						}

	const T operator()() const			{ return static_cast<T>(mersenne.Random());		}

private:
	mutable CRandomMersenne mersenne;

#endif
};

};	// namespace basis

#endif	
