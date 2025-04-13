#include <stdafx.h>
#include <string>
#include <Illumination/IlluminationHeader/Shape/IsoSurface_SRLGrid.h>
#include "SRL_Experiment.h"

using namespace std;
using namespace arithmetic;
using namespace fluidspace;
using namespace illuminationspace;

typedef illuminationspace::real T;

void SRL_Simple_VolumeRendering(const string filename, int frame) {
	int maxDepth = 6;

	illuminationspace::Scene* scene = new illuminationspace::Scene();
	XmlLibrary::importXmlFile(scene, "../SRL_Test/scene_SRL_XML.bdml");	

	vector<Primitive> vecPrims;

	InterfaceTracking_SRL3D srl;
	srl.load_compressed(fluidspace::get_concat_filename(filename, frame));

	bool isbdrySurface[3][2];
	for(int i=0;i<3;++i) for(int j=0;j<3;++j) isbdrySurface[i][j] = true;

	isbdrySurface[2][0] = true; 

	illuminationspace::Transform transform;

	BBox bbox(srl.srl->phi.getAABB());
	Vector3<T> center = bbox.center(); Vector3<T> extent = bbox.extent();

	Vector3<T> origin2center = bbox.min() + T(0.5f) * extent;
//	transform = arithmetic::translate(origin2center) * arithmetic::rotate(T(-10.0f), Vector3<T>(0,1,0)) * arithmetic::translate(-origin2center);

	IsoSurface_SRLGrid implicit(transform, *srl.srl, isbdrySurface);

	Spectrum _ex = Spectrum(5, 10, 160) / T(255);
	Water waterMat = Water(ConstantTexture<Spectrum>(Spectrum(1.0)), ConstantTexture<Spectrum>(Spectrum(1.0)), ConstantTexture<T>(1.33)
		, ConstantTexture<Spectrum>(log(_ex)), ConstantTexture<Spectrum>(log(Spectrum(1,1,1))));

	GeometricPrimitive p(implicit, waterMat);	vecPrims.push_back(p);	

	Shape sphere = illuminationspace::Sphere(transform, Vector3<T>(0.5*srl.get_Lx(),0.42*srl.get_Ly(),0.5*srl.get_Lz()), 0.1*srl.get_Ly());

	_ex = Spectrum(255, 0, 0) / T(255);
	Water sphereMat = Water(ConstantTexture<Spectrum>(Spectrum(1.0)), ConstantTexture<Spectrum>(Spectrum(1.0)), ConstantTexture<T>(1.33)
		, ConstantTexture<Spectrum>(log(_ex)), ConstantTexture<Spectrum>(log(Spectrum(1,1,1))));

	GeometricPrimitive sph(sphere, sphereMat);	vecPrims.push_back(sph);	

	scene->setAggregate(SimpleAggregate(vecPrims));
	scene->setBackground(Spectrum::white() * 0.8f);

	DirectLightingIntegrator* itg = new DirectLightingIntegrator(maxDepth);	
	scene->integrator = itg;

	dynamic_cast<ImageFilm*>(scene->camera->film)->filename = fluidspace::get_concat_filename("../SRL_Test/Image/srl", frame) + ".tga";

	scene->render();

	delete scene;
}

void SRL_Simple_VolumeRendering() {

	for(int i=0;i<200;++i) {
		SRL_Simple_VolumeRendering("../SRL_Test/BinaryData/Liquid3D _ InterfaceTracking_SRL3D", i);	
	}
	
}