#include <stdafx.h>
#include "SRL_Experiment.h"

using namespace std;
using namespace fluidspace;

void SRL_Experiment(Scene* scene) {
	//SRL_Simple_VolumeRendering();

	//SRL_FluidSimulation_DamBreaking(scene);
	
	//SRL_FluidSimulation_ObjectInteraction2D(scene);
	//SRL_FluidSimulation_ObjectInteraction(scene);
	
	//SRL_FluidSimulation_SmallScale2D(scene);
	//SRL_FluidSimulation_SmallScale(scene);	

	//SRL_FluidSimulation_3D_Demo_01(scene);		// falling bunny
	//SRL_FluidSimulation_3D_Demo_02(scene);		// falling armadillo
	//SRL_FluidSimulation_3D_Demo_03(scene);		// dam breaking
	//SRL_FluidSimulation_3D_Demo_04(scene);		// falling star
	//SRL_FluidSimulation_3D_Demo_05(scene);			// thin sheets

	//SRL_Binary2ImageData(scene,220);
	//SRL_Binary2OBJ(scene, 0, 600);

	SRL_Compare_Reinitialization_RotatingBunny_Infrequent(scene);	
	//SRL_Compare_Reinitialization_RotatingBunny(scene);	
	//SRL_Compare_Reinitialization_RotatingZalesak(scene);	
	//SRL_Compare_Reinitialization_2DTest(scene);
	
	//SRL_Compare_100x100_ZalesakTest(scene);
	//SRL_Compare_AdvectionTest_Rotation(scene);	

	//SRL_Measure_Scalability(scene);

#ifdef USE_OPENMP	
	int maxthreads = omp_get_max_threads();
	int numThreads = maxthreads;

	omp_set_num_threads(numThreads);
	log_info << endl << "OpenMP " << numThreads << " Threads are initiated" << endl << endl;
#endif

	//GridDistanceField3D<real> grid_field; 
	//grid_field.set("../Data/bunny.obj", 200, 200, 200, "../Data/bunny_phi");
	//grid_field.set("../Data/Armadillo.ply", 200, 200, 200, "../Data/Armadillo_phi", true, arithmetic::rotateY(180.0));	
	//grid_field.set("../Data/star3D.obj", 200, 200, 200, "../Data/star3D_phi");
	//grid_field.set("../Data/polyhedra_star.obj", 200, 200, 200, "../Data/polyhedra_star_phi");
}

void SRL_Binary2ImageData(Scene* scene, unsigned int frame) {
	InterfaceTracking_Viewer* viewer = new InterfaceTracking_Viewer();	
	viewer->enableSaveImage();
	viewer->load("../SRL_Test/Binary/Liquid3D _ InterfaceTracking_SRL3D");
	//viewer->load("../Output/FluidBinary _ WENO/RotatingBunnyWENO _ InterfaceTracking_SRL3D");

	if(frame != 0) {
		viewer->load_frame(frame);
	}

	scene->push_SceneObject(viewer);
}

void SRL_Binary2OBJ(Scene* scene, unsigned int frame_s = 0, unsigned int frame_e = 5000) {
	const string filename = output_folder + "SRLOBJ/srl_obj";

	InterfaceTracking_Viewer* viewer = new InterfaceTracking_Viewer();		
	viewer->load(output_folder + "FluidBinary/Liquid3D _ InterfaceTracking_SRL3D");

	for(unsigned int frame = frame_s; frame < frame_e ; frame++) {
		if(viewer->load_frame(frame) == false)
			continue;

		std::cout << "Load Frame " << frame << " is Done " << std::endl;

		char str[512];
		sprintf_s(str, std::string("_%04d").c_str(), frame);

		string obj_name = filename + std::string(str) + ".obj";

		viewer->exportOBJ(obj_name.c_str()); 
	}

//	cout << viewer->frame << " " << jump_frame << endl;
//	cout << frame_e << endl;
}
