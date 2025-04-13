#include <stdafx.h>
#include <Illumination/IlluminationHeader/Scene/Scene.h>

using namespace std;
using namespace basis;
using namespace illuminationspace;

Scene::~Scene() { 
	if(camera) delete camera; 
	if(imageSampler) delete imageSampler; 
	if(integrator) delete integrator; 
}

void Scene::init() {
	camera = NULL;
	imageSampler = NULL;
	integrator = NULL;	

	lights.clear();
	background = Spectrum::black();
	
	bbox.init();	
}

void Scene::render_serial() {
	Timer timer;

	if(integrator) {		
		integrator->preProcess(this);		// no problem
		timer.printTime("Time for Pre-Process");
		timer.reset();
	}

	Ray r;
	Spectrum radiance(0.0f);	
	ImageSample imageSample;	

	while(imageSampler->getNextSample(&imageSample)) {		// no problem
		integrate(r, radiance, imageSample);
		camera->film->addSample(imageSample, r, radiance, 1.0f);	// no problem
	}

	timer.printTime("Rendering : ");

	if(integrator) 
		integrator->postProcess(this);

	camera->film->writeImage();
	camera->film->print_debug_info();
}

void Scene::render_openmp() {
	Timer timer;

	if(integrator) {		
		integrator->preProcess(this);		// no problem
		timer.printTime("Time for Pre-Process");
		timer.reset();
	}

	const int nx = imageSampler->xPixelEnd-imageSampler->xPixelStart;
	const int ny = imageSampler->yPixelEnd-imageSampler->yPixelStart;
	const int n_pixels = nx*ny;

	Array2<Array<Ray>>			rays(nx,ny);
	Array2<Array<Spectrum>>		radiances(nx,ny);
	Array2<Array<ImageSample>>	samples(nx,ny);

	for(int i=0;i<nx;++i) { for(int j=0;j<ny;++j) {
		rays(i,j).resize(imageSampler->nSamplesPerPixel);
		radiances(i,j).resize(imageSampler->nSamplesPerPixel);
		samples(i,j).resize(imageSampler->nSamplesPerPixel);
	}}

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
	for(int i=0;i<nx;++i) { 
		for(int j=0;j<ny;++j) {
		imageSampler->getSamples(samples(i,j), i + imageSampler->xPixelStart, j + imageSampler->yPixelStart);

		for(size_t k=0;k<samples[i].size();++k)
			integrate(rays(i,j)(k), radiances(i,j)(k), samples(i,j)(k));
	}}

	if(integrator) 
		integrator->postProcess(this);

	timer.printTime("Rendering : "); timer.reset();

	for(size_t i=0;i<samples.size();++i)
		for(size_t k=0;k<samples[i].size();++k)
			camera->film->addSample(samples[i][k], rays[i][k], radiances[i][k], 1.0f);	// no problem

	camera->film->writeImage();
	camera->film->print_debug_info();

	timer.printTime("Write Film : "); timer.reset();	
}

void Scene::render()	{
	render_openmp();
	return;

#ifdef USE_OPENMP
	render_openmp();
#else	
	render_serial();
#endif	
}

void Scene::integrate(Ray& r, Spectrum& radiance, const ImageSample& imageSample) {
	camera->generateRay(imageSample, &r);						// no problem
	r.setTmin(Ray::Epsilon()); r.setTmax(Ray::Large());			// no problem

	integrator->getNextSamples();								// big problem _ sampler parallelization			
	radiance = Li(r);											// no problem

	radiance.clamp(0, LARGE_i);									// no problem	
	
	BSDF::FreeAll();											// now be cured
}

