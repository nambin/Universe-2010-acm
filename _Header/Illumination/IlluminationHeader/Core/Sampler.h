#ifndef _ILLUMINATION_SAMPLER_H_
#define _ILLUMINATION_SAMPLER_H_

#include <vector>
#include "../../../Basis/Basis/Environments.h"
#include "../../../Basis/XmlLibrary/XmlLibrary.h"
#include "../../IlluminationSpace.h"
#include "../Illumination/Random.h"
#include "../Illumination/Constants.h"

namespace illuminationspace {

using std::vector;

class Scene;

//////////////////////////////////////////////////////////////////////////
// Sampling Library
class ILLUMINATIONDLL SamplerLibrary {
public : 
	static void StratifiedSample1D(real *samp, int nSamples, bool jitter);
	static void StratifiedSample2D(real *samp, int nX, int nY, bool jitter);
	static void Shuffle(real *samp, int count, int dims);
	static void LatinHypercube(real *samples, int nSamples, int nDim);
};

class Sampler;

//////////////////////////////////////////////////////////////////////////
// Sampling
class ILLUMINATIONDLL Sampler_SingleThread {
public : 
	friend class Sampler;

protected : 
	enum SAMPLING_TYPE { RANDOM = 0, STRATIFIED_NOTJITTER, STRATIFIED_JITTER, LATINHYPERCUBE };

	Sampler_SingleThread();
	~Sampler_SingleThread();

	void allocate();

	void getNextSample();
	const unsigned int add1D(unsigned int num, SAMPLING_TYPE type);
	const unsigned int add2D(unsigned int num, unsigned int num_y, SAMPLING_TYPE type);

	void releaseMemory();
	void release();

	vector<unsigned int>	n1D, n2D, nX, nY;		// needed sample counting
	vector<SAMPLING_TYPE>	type1D, type2D;

	real					**oneD, **twoD;			// sampled value
};

class ILLUMINATIONDLL Sampler {
public : 
	typedef Sampler_SingleThread::SAMPLING_TYPE SAMPLING_TYPE;

	static const SAMPLING_TYPE RANDOM				= Sampler_SingleThread::RANDOM;
	static const SAMPLING_TYPE STRATIFIED_NOTJITTER = Sampler_SingleThread::STRATIFIED_NOTJITTER;
	static const SAMPLING_TYPE STRATIFIED_JITTER	= Sampler_SingleThread::STRATIFIED_JITTER;
	static const SAMPLING_TYPE LATINHYPERCUBE		= Sampler_SingleThread::LATINHYPERCUBE;

	Sampler() {}

#ifdef USE_OPENMP

	// PreProcess
	void allocate() { 
		const int max_thread = omp_get_max_threads();
		for(int k=0;k<max_thread;++k)
			sampler[k].allocate(); 
	}
	
	const unsigned int add1D(unsigned int num, SAMPLING_TYPE type = RANDOM)	{ 
		const int max_thread = omp_get_max_threads();
		unsigned int ret;
		for(int k=0;k<max_thread;++k)
			ret = sampler[k].add1D(num,type);			
		return ret;
	}
	const unsigned int add2D(unsigned int num, unsigned int num_y = 0, SAMPLING_TYPE type = RANDOM) { 
		const int max_thread = omp_get_max_threads();
		unsigned int ret;
		for(int k=0;k<max_thread;++k)
			ret = sampler[k].add2D(num,num_y,type);	
		return ret;
	}

	// Run-Time Process
	void getNextSample() { sampler[omp_get_thread_num()].getNextSample(); }

	const real get1D(size_t i, size_t j) const { return sampler[omp_get_thread_num()].oneD[i][j]; }
	const real get2D(size_t i, size_t j) const { return sampler[omp_get_thread_num()].twoD[i][j]; }

	const unsigned int num1D(size_t i) const { return sampler[omp_get_thread_num()].n1D[i]; }
	const unsigned int num2D(size_t i) const { return sampler[omp_get_thread_num()].n2D[i]; }

protected : 
	Sampler_SingleThread	sampler[NUM_OPENMP_THREADS_MAX];

#else 

	// PreProcess
	void allocate() { sampler.allocate(); }
	
	const unsigned int add1D(unsigned int num, SAMPLING_TYPE type = RANDOM)							{ return sampler.add1D(num,type);		}
	const unsigned int add2D(unsigned int num, unsigned int num_y = 0, SAMPLING_TYPE type = RANDOM) { return sampler.add2D(num,num_y,type); }

	// Run-Time Process
	void getNextSample()																			{ sampler.getNextSample();				}

	const real get1D(size_t i, size_t j) const { return sampler.oneD[i][j]; }
	const real get2D(size_t i, size_t j) const { return sampler.twoD[i][j]; }

	const unsigned int num1D(size_t i) const { return sampler.n1D[i]; }
	const unsigned int num2D(size_t i) const { return sampler.n2D[i]; }

protected : 
	Sampler_SingleThread	sampler;

#endif	
};

class ILLUMINATIONDLL Sampler_Offset {
public : 
	typedef Array<unsigned int>::index_type index_type;

	Sampler_Offset() {}

#ifdef USE_OPENMP

	// PreProcess
	void resize(index_type num) { 
		const int max_thread = omp_get_max_threads();
		for(int k=0;k<max_thread;++k)
			offset[k].resize(num);		
	}

	void set(index_type i, unsigned int val) { 
		const int max_thread = omp_get_max_threads();
		for(int k=0;k<max_thread;++k)
			offset[k][i] = val;
	}
	unsigned int& operator[](index_type i) { return offset[omp_get_thread_num()][i]; }

	// Run-Time Process	
	const unsigned int	operator()(index_type i) const	{ return offset[omp_get_thread_num()][i]; }	
	const unsigned int	operator[](index_type i) const	{ return offset[omp_get_thread_num()][i]; }

private :
	Array<unsigned int> offset[NUM_OPENMP_THREADS_MAX];

#else

	// PreProcess
	void resize(index_type num)					{ offset.resize(num);	}
	void set(index_type i, unsigned int val)	{ offset[i] = val;		}

	// Run-Time Process
	const unsigned int	operator()(index_type i) const	{ return offset[i]; }
	const unsigned int	operator[](index_type i) const	{ return offset[i]; }

private :
	Array<unsigned int> offset;

#endif
};

//template<typename T> 
//void Sampler::requestSamples(T* ptr, const Scene* scene) {
//	if(ptr != NULL) 
//		ptr->requestSamples(this, scene);
//	allocate();	
//}

};		// illuminationspace

#endif