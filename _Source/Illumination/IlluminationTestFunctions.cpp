#include <stdafx.h>
#include <Illumination/IlluminationHeader/IlluminationTestFunctions.h>
#include <Principle/Geometry/Levelset/Levelset3D.h>
#include <Principle/Geometry/SRL/SRL3D.h>
#include <Principle/Geometry/DistanceFieldLibrary/DistanceFields.h>
#include <Illumination/IlluminationHeader/Shape/IsoSurface_SRLGrid.h>

using namespace std;
using namespace basis;
using namespace geometryspace;

namespace illuminationspace {

Scene *scene = NULL;

const real	lightRadiance = 20.0f;
const int	sqrtLightSample = 1;
const int	maxDepth = 6;

const real	L = 2.0f;
const real	Lx = L * 0.5f, Ly = L * 0.5f, Lz = L * 0.5f;
const Vector3<real> trans(Lx, Ly, 2*Lz);

void test() {	
	int numThreads = 8;

#ifdef USE_OPENMP		
	omp_set_num_threads(omp_get_max_threads());
	log_info << endl << "OpenMP " << omp_get_max_threads() << " Threads are initiated" << endl << endl;
#endif
	
//	makeScene_CornellBox();
//	makeScene_Levelset_ImplicitRendering();
	makeScene_SRL_ImplicitRendering();

	scene->render();
}

void makeScene_Levelset_ImplicitRendering() {
	vector<Primitive> vecPrims;
	set_Levelset_ImplicitPrimitive(vecPrims);	

	scene = new Scene();
	XmlLibrary::importXmlFile(scene, "../Illumination/Data/Scene/SceneXML.bdml");	

	scene->setBackground(Spectrum::white() * 0.8f);
	scene->setAggregate(SimpleAggregate(vecPrims));
	scene->integrator = new DirectLightingIntegrator(maxDepth);	
}

void makeScene_SRL_ImplicitRendering() {
	scene = new Scene();
	XmlLibrary::importXmlFile(scene, "../Illumination/Data/Scene/SceneXML.bdml");	

	vector<Primitive> vecPrims;
	set_SRL_ImplicitPrimitive(vecPrims);
	set_CornellBox(vecPrims);

	vector<Light> lts;
	set_Lights(vecPrims, lts); scene->setLights(lts);
	scene->setAggregate(SimpleAggregate(vecPrims));
	scene->setBackground(Spectrum::black());
	
	//DirectLightingIntegrator* itg = new DirectLightingIntegrator(maxDepth);	
	
	PhotonMapIntegrator* itg = new PhotonMapIntegrator(maxDepth);
	itg->init_UnbiasedPhotonEmitter(1.0e+6);
	itg->init_BasicPhotonGather(1000, 0.05f * L);
	itg->directPhotonRendering = false;

	scene->integrator = itg;

	itg->set_lightIntegrator(new LightIntegrator_UniformSampleAll());	
}

void set_Levelset_ImplicitPrimitive(std::vector<Primitive> &vecPrims) {	
	int dim[3] = { 100, 100, 100 };

	BBox bbox; set_BBox(bbox);
	Vector3<real> center = bbox.center(); Vector3<real> extent = bbox.extent();
	
	Vector3<real> _center(0.3*L, 0.6*L, 0.3*L);	
	real radius = 0.25 * L;

	Levelset3D<real> lv; lv.set(dim[0]-1,dim[1]-1,dim[2]-1,extent[0]/(dim[0]-1));
	lv.setPhi(WaterDropDistanceField3D<real>(_center, radius, 0.2*L));
	
	bool isbdrySurface[3][2];
	for(int i=0;i<3;++i) for(int j=0;j<3;++j) isbdrySurface[i][j] = true;
	
	isbdrySurface[2][0] = true; 

	Vector3<real> origin2center = bbox.min() + real(0.5f) * extent;

	Transform identity;
	identity = arithmetic::translate(origin2center) * arithmetic::rotate(real(-10.0f), Vector3<real>(0,1,0)) * arithmetic::translate(-origin2center);
	
	IsoSurface_RegularGrid implicit(identity, lv.phi, isbdrySurface);

	//Spectrum _ex = Spectrum(5, 10, 160) / real(255);
	Spectrum _ex = Spectrum(255, 255, 255) / real(255);
	Water waterMat = Water(ConstantTexture<Spectrum>(Spectrum(1.0)), ConstantTexture<Spectrum>(Spectrum(1.0)), ConstantTexture<real>(1.10)
						 , ConstantTexture<Spectrum>(log(_ex)), ConstantTexture<Spectrum>(log(Spectrum(1,1,1))));

	GeometricPrimitive p(implicit, waterMat);	vecPrims.push_back(p);	
}

void set_SRL_ImplicitPrimitive(std::vector<Primitive> &vecPrims) {	
	int dim[3] = { 100, 100, 100 };

	BBox bbox; set_BBox(bbox);
	Vector3<real> center = bbox.center(); Vector3<real> extent = bbox.extent();

	Vector3<real> _center(0.3*L, 0.6*L, 0.3*L);	
	real radius = 0.25 * L;

	SRL3D<real> srl; srl.set(dim[0]-1,dim[1]-1,dim[2]-1,extent[0]/(dim[0]-1));
	srl.setPhi(WaterDropDistanceField3D<real>(_center, radius, 0.2*L));

	bool isbdrySurface[3][2];
	for(int i=0;i<3;++i) for(int j=0;j<3;++j) isbdrySurface[i][j] = true;

	isbdrySurface[2][0] = true; 

	Transform transform;
	Vector3<real> origin2center = bbox.min() + real(0.5f) * extent;
//	transform = arithmetic::translate(origin2center) * arithmetic::rotate(real(-10.0f), Vector3<real>(0,1,0)) * arithmetic::translate(-origin2center);

	IsoSurface_SRLGrid implicit(transform, srl, isbdrySurface);

	//Spectrum _ex = Spectrum(5, 10, 160) / real(255);
	Spectrum _ex = Spectrum(150, 150, 160) / real(255);
	Water waterMat = Water(ConstantTexture<Spectrum>(Spectrum(1.0)), ConstantTexture<Spectrum>(Spectrum(1.0)), ConstantTexture<real>(1.33)
						 , ConstantTexture<Spectrum>(log(_ex)), ConstantTexture<Spectrum>(log(Spectrum(1,1,1))));

	GeometricPrimitive p(implicit, waterMat);	vecPrims.push_back(p);	
}

void makeScene_CornellBox() {
	vector<Primitive> vecPrims;
	vector<Light> lts;

	scene = new Scene();		

	XmlLibrary::importXmlFile(scene, "../Illumination/Data/Scene/SceneXML.bdml");	

	set_CornellBox(vecPrims); set_CornellSphere(vecPrims);
	set_Lights(vecPrims, lts); scene->setLights(lts);	

	//scene->setAggregate(KdTreeAccel(vecPrims));
	//scene->setAggregate(GridAccel(vecPrims));
	scene->setAggregate(SimpleAggregate(vecPrims));

	DirectLightingIntegrator* itg = new DirectLightingIntegrator(maxDepth);
	//PhotonMapIntegrator* itg = new PhotonMapIntegrator(maxDepth);
	//itg->init_UnbiasedPhotonEmitter(1.0e+6);
	//itg->init_BasicPhotonGather(1000, 0.05f * L);
	//itg->directPhotonRendering = false;

	itg->set_lightIntegrator(new LightIntegrator_UniformSampleAll());
	//itg->set_lightIntegrator(new LightIntegrator_UniformSampleOne());
	scene->integrator = itg;
}

void set_CornellBox(std::vector<Primitive> &vecPrims) {
	Transform identity;

	Shape sQuadBack = Parallelogram(identity, trans+Vector3<real>(-Lx, -Ly, -2.0*Lz), Vector3<real>(2.0*Lx, 0.0, 0.0), Vector3<real>(0.0, 2.0*Ly, 0.0));
	Shape sQuadSide1 = Parallelogram(identity, trans+Vector3<real>(-Lx, -Ly, 0.0f), Vector3<real>(0.0, 0.0, -2.0*Lz), Vector3<real>(0.0, 2.0*Ly, 0.0));
	Shape sQuadSide2 = Parallelogram(identity, trans+Vector3<real>( Lx, -Ly, -2.0f*Lz), Vector3<real>(0.0, 0.0, 2.0*Lz), Vector3<real>(0.0, 2.0*Ly, 0.0));
	Shape sQuadUp = Parallelogram(identity, trans+Vector3<real>(-Lx,  Ly, -2.0f*Lz), Vector3<real>(2.0*Lx, 0.0, 0.0), Vector3<real>(0.0, 0.0, 2.0*Lz));
	Shape sQuadBottom = Parallelogram(identity, trans+Vector3<real>(-Lx, -Ly, 0.0f), Vector3<real>(2.0*Lx, 0.0, 0.0), Vector3<real>(0.0, 0.0, -2.0*Lz));

	Material mQuadBack = Matte(ConstantTexture<Spectrum>(Spectrum(0.7, 0.7, 0.7) * 0.75), ConstantTexture<real>(0.0));
	Material mQuadSide1 = Matte(ConstantTexture<Spectrum>(Spectrum(0.4, 0.2, 0.2) * 2), ConstantTexture<real>(0.0));
	Material mQuadSide2 = Matte(ConstantTexture<Spectrum>(Spectrum(0.2, 0.2, 0.4) * 2), ConstantTexture<real>(0.0));
	Material mQuadUp = Matte(ConstantTexture<Spectrum>(Spectrum(0.7, 0.7, 0.7) * 0.75), ConstantTexture<real>(0.0));
	Material mQuadBottom = Matte(ConstantTexture<Spectrum>(Spectrum(0.7, 0.7, 0.7) * 0.75), ConstantTexture<real>(0.0));

	GeometricPrimitive pQuadBack = GeometricPrimitive(sQuadBack, mQuadBack);
	GeometricPrimitive pQuadSide1 = GeometricPrimitive(sQuadSide1, mQuadSide1);
	GeometricPrimitive pQuadSide2 = GeometricPrimitive(sQuadSide2, mQuadSide2);
	GeometricPrimitive pQuadUp = GeometricPrimitive(sQuadUp, mQuadUp);
	GeometricPrimitive pQuadBottom = GeometricPrimitive(sQuadBottom, mQuadBottom);
	
	vecPrims.push_back(pQuadBack);
	vecPrims.push_back(pQuadSide1);
	vecPrims.push_back(pQuadSide2);
	vecPrims.push_back(pQuadBottom);
	//vecPrims.push_back(pQuadUp);
}

void set_CornellSphere(std::vector<Primitive> &vecPrims) {
	Transform identity;

	Shape sSphere1 = Sphere(identity, trans+Vector3<real>(-0.45*Lx,-Ly+0.36*Ly,-1.0*Lz), 0.36*Ly);
	Shape sSphere2 = Sphere(identity, trans+Vector3<real>( 0.47*Lx,-Ly+0.36*Ly,-0.8*Lz), 0.36*Ly);

	Material mSphere1 = Mirror(ConstantTexture<Spectrum>(Spectrum(0.8)));
	Material mSphere2 = Glass(ConstantTexture<Spectrum>(Spectrum(1.0)), ConstantTexture<Spectrum>(Spectrum(1.0)), ConstantTexture<real>(2.3));	

	GeometricPrimitive pSphere1 = GeometricPrimitive(sSphere1, mSphere1);
	GeometricPrimitive pSphere2 = GeometricPrimitive(sSphere2, mSphere2);

	vecPrims.push_back(pSphere1); vecPrims.push_back(pSphere2);		
}

void set_Lights(std::vector<Primitive> &vecPrims, std::vector<Light> &lts) {
	////////////////////////////////////////////////////////////////////////////////////////
	// Set Lights
	real ratioX = 0.45f;
	real ratioZ = 0.30f;
	real centerZ = -1.2f*Lz;

	Sampler::SAMPLING_TYPE		type = Sampler::STRATIFIED_NOTJITTER;
	if(sqrtLightSample != 1)	type = Sampler::STRATIFIED_JITTER;

	Shape sLight = Parallelogram(Transform(), trans+Vector3<real>(-ratioX*Lx, 0.990f*Ly, centerZ-ratioZ*Lz), Vector3<real>(2.0f*ratioX*Lx, 0.0, 0.0), Vector3<real>(0.0, 0.0, 2.0*ratioZ*Lz));
	AreaLight areaLight = AreaLight(Transform(), Spectrum(lightRadiance), sLight, sqrtLightSample, sqrtLightSample, type);
	//AreaLight areaLight = AreaLight(Transform(), Spectrum::blue(), sLight, sqrtLightSample, sqrtLightSample, type);
	GeometricPrimitive pLight = GeometricPrimitive(sLight, BlackMaterial(), areaLight);	

	vecPrims.push_back(pLight);	
	lts.push_back(areaLight);
}

void set_BBox(BBox &bbox) {
	Vector3<real> _min(-Lx, -Ly, -2.0f*Lz), _max(Lx, Ly, 0.0f);
	bbox.set(trans+_min, trans+_max);
}

};