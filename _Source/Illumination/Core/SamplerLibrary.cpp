#include <stdafx.h>
#include <Illumination/IlluminationHeader/Core/Sampler.h>

using namespace illuminationspace;

void SamplerLibrary::StratifiedSample1D(real *samp, int nSamples, bool jitter) {
	real invTot = 1.0f / static_cast<real>(nSamples);
	for(int i=0;i<nSamples;++i) {
		real j = jitter ? random() : 0.5f;
		*samp++ = (i+j) * invTot;
	}
}

void SamplerLibrary::StratifiedSample2D(real *samp, int nX, int nY, bool jitter) {
	real dx = 1.0f/static_cast<real>(nX), dy = 1.0f/static_cast<real>(nY);
	for(int i=0;i<nX;++i) {
		for(int j=0;j<nY;++j) {
			real jx = jitter ? random() : 0.5f;
			real jy = jitter ? random() : 0.5f;
			*samp++ = (i+jx) * dx;
			*samp++ = (j+jy) * dy;
		}
	}
}

void SamplerLibrary::Shuffle(real *samp, int count, int dims) {
	for(int i=0;i<count;++i) {
		int other = ((unsigned int)random.getInt() % count);
		for(int j=0;j<dims;++j)
			std::swap(samp[dims*i+j], samp[dims*other+j]);
	}
}

void SamplerLibrary::LatinHypercube(real *samples, int nSamples, int nDim) {
	real delta = 1.0f / static_cast<real>(nSamples);
	for(int i=0;i<nSamples;++i) {
		for(int j=0;j<nDim;++j) {
			samples[nDim*i+j] = (i+random()) * delta;
		}
	}

	for(int i=0;i<nDim;++i) {
		for(int j=0;j<nSamples;++j) {
			int other = ((unsigned int)random.getInt() % nSamples);
			std::swap(samples[nDim*j+i], samples[nDim*other+i]);
		}
	}
}
