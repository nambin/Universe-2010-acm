#include <stdafx.h>
#include <Basis/MemoryManagement/MemoryManagement.h>
#include <Illumination/IlluminationHeader/Core/Sampler.h>
#include <Illumination/IlluminationHeader/Integrator/SurfaceIntegrator.h>
#include <Illumination/IlluminationHeader/Integrator/VolumeIntegrator.h>

using namespace illuminationspace;

Sampler_SingleThread::Sampler_SingleThread() { oneD = twoD = NULL; }
Sampler_SingleThread::~Sampler_SingleThread() { release(); }

void Sampler_SingleThread::allocate() {
	releaseMemory();

	int nPtrs = static_cast<int>(n1D.size() + n2D.size());
	if(!nPtrs) {
		oneD = twoD = NULL;
		return;
	}
	oneD = (real**)basis::AllocAligned(nPtrs * sizeof(real*));
	twoD = oneD + n1D.size();

	int totalSamples = 0;
	for(unsigned int i=0;i<n1D.size();++i) totalSamples += n1D[i];
	for(unsigned int i=0;i<n2D.size();++i) totalSamples += 2 * n2D[i];

	real *mem = (real*)basis::AllocAligned(totalSamples * sizeof(real));
	for(unsigned int i=0;i<n1D.size();++i) {
		oneD[i] = mem;
		mem += n1D[i];
	}
	for(unsigned int i=0;i<n2D.size();++i) {
		twoD[i] = mem;
		mem += 2 * n2D[i];
	}
}

void Sampler_SingleThread::releaseMemory() {
	if(oneD == NULL && twoD == NULL)
		return;

	int totalSamples = 0;
	for(unsigned int i=0;i<n1D.size();++i) totalSamples += n1D[i];
	for(unsigned int i=0;i<n2D.size();++i) totalSamples += 2 * n2D[i];
	if(totalSamples != 0)
		basis::FreeAligned(oneD[0]);

	int nPtrs = static_cast<int>(n1D.size() + n2D.size());
	if(nPtrs != 0)
		basis::FreeAligned(oneD);

	oneD = twoD = NULL;
}

void Sampler_SingleThread::release() {	
	releaseMemory();
	n1D.clear();	n2D.clear();
	nX.clear();		nY.clear();
	type1D.clear();	type2D.clear();	
}

void Sampler_SingleThread::getNextSample() {
	for(unsigned int i=0;i<n1D.size();++i)	{
		switch(type1D[i]) {
			case RANDOM :
				for(unsigned int j=0;j<n1D[i];++j)
					oneD[i][j] = random();
				break;

			case STRATIFIED_NOTJITTER :
				SamplerLibrary::StratifiedSample1D(oneD[i], n1D[i], false);
				break;

			case STRATIFIED_JITTER :
				SamplerLibrary::StratifiedSample1D(oneD[i], n1D[i], true);
				break;

			case LATINHYPERCUBE :
				SamplerLibrary::LatinHypercube(oneD[i], n1D[i], 1);
				break;
		}	
	}

	for(unsigned int i=0;i<n2D.size();++i) {
		switch(type2D[i]) {
			case RANDOM :
				for(unsigned int j=0;j<n2D[i];++j) {
					twoD[i][2*j  ] = random();
					twoD[i][2*j+1] = random();
				}
				break;

			case STRATIFIED_NOTJITTER :
				SamplerLibrary::StratifiedSample2D(twoD[i], nX[i], nY[i], false);
				break;

			case STRATIFIED_JITTER :
				SamplerLibrary::StratifiedSample2D(twoD[i], nX[i], nY[i], true);
				break;

			case LATINHYPERCUBE : 
				SamplerLibrary::LatinHypercube(twoD[i], n2D[i], 2);
				break;
		}
	}
}

const unsigned int Sampler_SingleThread::add1D(unsigned int num, SAMPLING_TYPE type) {
	n1D.push_back(num);
	type1D.push_back(type);
	return static_cast<unsigned int>(n1D.size()-1);
}

const unsigned int Sampler_SingleThread::add2D(unsigned int num, unsigned int num_y, SAMPLING_TYPE type) {

	if(type == STRATIFIED_JITTER || type == STRATIFIED_NOTJITTER) {
		n2D.push_back(num*num_y);
		nX.push_back(num);
		nY.push_back(num_y);
		type2D.push_back(type);
	}
	else {
		n2D.push_back(num);
		nX.push_back(num);
		nY.push_back(0);
		type2D.push_back(type);
	}

	return static_cast<unsigned int>(n2D.size()-1);
}