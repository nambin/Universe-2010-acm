#ifndef _SCENE_NODE_ADDER_H_
#define _SCENE_NODE_ADDER_H_

using namespace controlspace;
using namespace arithmetic;
using namespace geometryspace;
using namespace shaderspace;
using namespace nature;
using namespace fluidspace;
using namespace rigidspace;
using namespace springspace;

int numThreads = 8;

namespace controlspace {

void InitLights(Scene* scene) {

	RGBA ambient = RGBA(1.0f, 1.0f, 1.0f, 1.0f) * 0.05f * 1;
	RGBA diffuse = RGBA(1.0f, 1.0f, 1.0f, 1.0f) * 0.05f * 1;
	RGBA specular = RGBA(1.0f, 1.0f, 1.0f, 1.0f);

	// Set Lights
	Light light0, light1, light2, light3;

	// Cloth Scene
	//light0.setAmbient(ambient);
	light0.setDiffuse(diffuse * 1.0f);
	light0.setSpecular(specular * 1.0f);
	light0.setPosition(Vector3f(-1.0f,5.0f,-1.5f) * 10.0f);
	scene->push_SceneLight(light0);

	//light1.setAmbient(ambient);
	//light1.setDiffuse(diffuse);
	//light1.setSpecular(specular);
	//light1.setPosition(Vector3f(0,5.0f,1.5f) * 0.1f);
	//scene->push_SceneLight(light1);

	//light2.setAmbient(ambient);
	//light2.setDiffuse(diffuse);
	//light2.setSpecular(specular);
	//light2.setPosition(Vector3f(0,1.0f,-5.0f) * 0.1f);
	//scene->push_SceneLight(light2);

	// Enright Test	
	light2.setDiffuse(diffuse * 5.5f);
	light2.setSpecular(specular * 0.1f);
	light2.setPosition(Vector3f(1.0f,-1.0f,-1.0f) * 10.0f);
	scene->push_SceneLight(light2);

	// Fluid
	light3.setAmbient(ambient * 10.0f);
	light3.setDiffuse(diffuse * 10.0f);
	light3.setSpecular(specular * 1.0f);
	light3.setPosition(Vector3f(8.0f,8.0f,8.0f) * 10.0f);
	scene->push_SceneLight(light3);

	scene->setLightSize(0.01f);
	scene->set_SceneLightVisible(true);
}

void AddOpenGLTest(Scene * scene) {
	TestOpenGLScene *test_scene = new TestOpenGLScene;
	scene->push_SceneObject(test_scene);
}

void AddFluidRelated(Scene* scene) {
	typedef fluidspace::real T;

	//WENO_Adaptive<double>::validate_weights_uniform_weno(0.3);

	//double h[6] = { 0.1, 0.15, 0.17, 0.10, 0.08, 0.13 };	//for(int i=0;i<6;++i) h[i] = 1;
	//WENO_Adaptive<double>::validate_weights_except_smoothness(h);

	//WENO_Adaptive<double>::validate_values_uniform_weno(0.131, Sin<double>(5, 1, 3), 0, 10 * PI, 10);

	WENO3rd_Adaptive<double>::validate_weights_uniform_weno(0.31);

	double h[4] = { 0.1, 0.15, 0.17, 0.12 };	//for(int i=0;i<6;++i) h[i] = 1;
	WENO3rd_Adaptive<double>::validate_weights_except_smoothness(h);

	WENO3rd_Adaptive<double>::validate_values_uniform_weno(0.131, Sin<double>(5, 1, 3), 0, 10 * PI, 30);

	exit(1);

	//double h;
	//int nX, nY, nZ;
	//int resolution = 32;
	//int endFrame = 200;
	//double dt = 1.0f/30.0f;	

//	FluidScene *fluidscene = new FluidScene(dt);
//	fluidscene->initFluidSceneIntegrator(new FluidSceneIntegrator_StaticRigidObject(dt));	
//	fluidscene->setEndFrame(endFrame);
////	fluidscene->enableSaveBinary(500);
////	fluidscene->enableSaveRenderObjectsOnly(500);
//	fluidscene->enableSaveImage();		
//	scene->push_SceneObject(fluidscene);
//
//	Liquid2D *liquid2D = new Liquid2D(resolution*3, resolution * 2, 1.0f/resolution, dt, 2.0f, 2.0f);
//	liquid2D->setPhi(Grid2r());
//	fluidscene->initFluidSimulator(liquid2D);

	//FluidScene_Objects2D* object2D = new FluidScene_Objects2D();
	//object2D->push_DynamicObject(StaticRigidBody2D<T>(Circle<T>(Vector2<T>(0.5f, 0.6f), 0.13f)));
	//fluidscene->initFluidScene_Objects(object2D);
	
	//nX = nZ = 64 * 2; nY = 64 * 2;
	//h = 0.05 / double(nX);

	//Liquid3D *liquid3D = new Liquid3D(nX, nY, nZ, h, dt, 3.5f, 2.0f);
	//liquid3D->setPhi(Grid3r());
	//fluidscene->initFluidSimulator(liquid3D);

	//FluidScene_Objects3D* object3D = new FluidScene_Objects3D();
	//object3D->push_DynamicObject(StaticRigidBody3D<T>(Sphere<T>(Vector3<T>(nX*h*0.5f, nY*h*0.42f, nZ*h*0.5f), nY*h*0.10f)));
	//fluidscene->initFluidScene_Objects(object3D);

	//Smoke2D *smoke2D = new Smoke2D(resolution, resolution * 3 / 2, 1.0f / resolution, 1.0f/30.0f, 2.0f);
	//smoke2D->setDensity(Grid2r());	
	//fluidscene->initFluidSimulator(smoke2D);

	//Smoke3D *smoke3D = new Smoke3D(resolution, resolution * 3 / 2, resolution, 1.0f/resolution, 1.0f/30.0f, 2.0f);
	//smoke3D->setDensity(Grid3r());
	//fluidscene->initFluidSimulator(smoke3D);

//	InterfaceTrackingTest *test = new InterfaceTrackingTest();
	//test->setExperiment_Translation_2D();
//	test->setExperiment_Zalesak_2D();	
	//test->setExperiment_Translation_3D();
	//test->setExperiment_Zalesak_3D();	
	//test->setExperiment_Enright_3D();	
	//test->enableSaveBinary();
	//test->enableSaveImage();
//	scene->push_SceneObject(test);

//	InterfaceTracking_Viewer* viewer = new InterfaceTracking_Viewer();	
//	viewer->enableSaveImage();
////	viewer->load(output_folder + "FluidBinary/Zalesak 2D _ IntefaceTracking_SRL2D");
////	viewer->load(output_folder + "FluidBinary _ SLC _ Zalesak _ 256/InterfaceTracking_SLC3D");
////	viewer->load(output_folder + "FluidBinary/Liquid2D _ IntefaceTracking_SRL2D");
//	viewer->load(output_folder + "FluidBinary/Liquid3D _ InterfaceTracking_SRL3D");
//	scene->push_SceneObject(viewer);

	//ReinitializationAccuracyTest* test = new ReinitializationAccuracyTest();
	//test->setExperiment_Russo();
	//scene->push_SceneObject(test);

#ifdef USE_OPENMP	
	int maxthreads = omp_get_max_threads();
	int numThreads = maxthreads;

	omp_set_num_threads(numThreads);
	log_info << endl << "OpenMP " << numThreads << " Threads are initiated" << endl << endl;
#endif
}

void AddMassSpringRelated(Scene *scene) {

	typedef springspace::real real;

	real len = 0.3f;
	real dt0 = 1.0f/500.0f;
	real cfl = 1.0f;
	size_t endframe = 360;

	real stiffness = 1.0e+6;
	real damping = 3.0e+2;

	//// Quad Mesh
	//QuadDynamicMesh<real>* mesh = new QuadDynamicMesh<real>(35, 35, len, len);
	//mesh->setMass(500);
	//mesh->make_StretchSpring(new SpringFunctor_ElasticDamping<real>(stiffness, damping));
	//mesh->make_ShearSpring(new SpringFunctor_ElasticDamping<real>(stiffness, damping));

	//mesh->setTwoPointGrabConstraint();

	// Tri Mesh
	TriangleManipulator triman;
	triman.setFlag_CDT(true);
	//triman.setArealimit(1.0e-6 * 20);
	triman.setArealimit(1.0e-6 * 50);
	triman.makeRectangleMesh(Vector2<double>(0), Vector2<double>(len, len));

	TriDynamicMesh<real> *mesh = new TriDynamicMesh<real>(triman.mesh);
	mesh->setMass(500);
	mesh->make_StretchSpring(new SpringFunctor_ElasticDamping<real>(stiffness, damping));
	mesh->fixVertex(0);	mesh->fixVertex(1);

	mesh->transformMesh(rotateX(real(-90)));
	mesh->transformMesh(translate(Vector3<real>(-len*0.5f, len*0.7f, len*0.5f)));
	mesh->transformMesh(rotateY(real(180)));

	// Cloth Integrator
	ClothSceneIntegrator* integrator = new ClothSceneIntegrator_Basic(dt0);
	integrator->init_ImplicitEuler_Integrator(*mesh, cfl);
	//integrator->init_Bridson03_Integrator(*mesh, cfl);
	integrator->init_Bridson03_ObjectCollision(dt0);	

	// Object Scene
	ClothScene_Objects* objects = new ClothScene_Objects();

	StaticRigidBody3D<real> sphere(Sphere<real>(Vector3<real>(0,len/6.0f,0), len/6.0f));
	sphere.setFrictionCoefficient(1.0f);
	objects->push_DynamicObject(sphere);

	StaticRigidBody3D<real> plane(Plane<real>(Vector3<real>(0,1,0), 0));
	plane.setFrictionCoefficient(1.0f);
	objects->push_DynamicObject(plane);	

	ClothScene* clothScene = new ClothScene(1.0f/60.0f);
	//ClothScene* clothScene = new ClothScene(dt0);
	scene->push_SceneObject(clothScene);
//	clothScene->enableSaveImage();
	clothScene->setEndFrame(endframe);	
	
	clothScene->initDynamicMesh(mesh);
	clothScene->initClothScene_Objects(objects);
	clothScene->initClothSceneIntegrator(integrator);
}

void AddGeometryRelated(Scene* scene) {
	//TriangleManipulator triman;
	//triman.setFlag_CDT(true);
	//triman.setArealimit(0.0005f);
	//Vector2<double> lb = Vector2<double>(-0.5f,-0.5f);
	//triman.makeRectangleMesh(lb, -lb);

	//TriangleMeshViewer *viewer = new TriangleMeshViewer(triman.mesh);
	//scene->push_SceneObject(viewer);

	//ObjectDistanceField3D<double> bunny_dist("../Data/OBJs/bunny.obj");

	MeshViewer *meshviewer = new MeshViewer();
	scene->push_SceneObject(meshviewer);	
//	meshviewer->load("../Data/bunny.obj");		
//	meshviewer->load("../Data/dragon_high.obj");		
//	meshviewer->load("../Data/venusm.obj");	
	meshviewer->load("../Data/Armadillo.ply", false);	

	//QuadTreeViewer *quadtreeviewer = new QuadTreeViewer();
	//scene->push_SceneObject(quadtreeviewer);

	//OctreeViewer *octreeviewer = new OctreeViewer();
	//scene->push_SceneObject(octreeviewer);

	//ShaderViewer *viewer = new ShaderViewer();
	//scene->push_SceneObject(viewer);

	//BVH<size_t, float> bvh;	
}

}

#endif