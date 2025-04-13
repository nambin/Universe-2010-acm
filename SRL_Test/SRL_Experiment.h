#pragma once

#pragma warning(disable:4251)

#ifdef _WIN32
#ifdef SRL_TEST
#define SRLDLL __declspec(dllexport)
#else
#define SRLDLL __declspec(dllimport)
#endif
#define DLLEXPORT __declspec(dllexport)
#else
#define SRLDLL
#define DLLEXPORT
#endif

SRLDLL void SRL_Experiment(Scene* scene);

SRLDLL void SRL_Binary2ImageData(Scene* scene, unsigned int frame = 0);
SRLDLL void SRL_Binary2OBJ(Scene* scene, unsigned int frame_s, unsigned int frame_e);

SRLDLL void SRL_Simple_VolumeRendering();

SRLDLL void SRL_FluidSimulation_DamBreaking(Scene*);

SRLDLL void SRL_FluidSimulation_ObjectInteraction(Scene* scene);
SRLDLL void SRL_FluidSimulation_ObjectInteraction2D(Scene* scene);

SRLDLL void SRL_FluidSimulation_SmallScale(Scene* scene);
SRLDLL void SRL_FluidSimulation_SmallScale2D(Scene* scene);

SRLDLL void SRL_FluidSimulation_3D_Demo_01(Scene*);		// bunny falling
SRLDLL void SRL_FluidSimulation_3D_Demo_02(Scene*);		// armadillo falling
SRLDLL void SRL_FluidSimulation_3D_Demo_03(Scene*);		// dam breaking
SRLDLL void SRL_FluidSimulation_3D_Demo_04(Scene*);		// star Falling
SRLDLL void SRL_FluidSimulation_3D_Demo_05(Scene*);		// thin sheets

SRLDLL void SRL_Compare_Reinitialization_RotatingBunny(Scene*);
SRLDLL void SRL_Compare_Reinitialization_RotatingBunny_Infrequent(Scene*);

SRLDLL void SRL_Compare_Reinitialization_RotatingZalesak(Scene*);
SRLDLL void SRL_Compare_Reinitialization_2DTest(Scene*);

SRLDLL void SRL_Compare_AdvectionTest_Rotation(Scene*);
SRLDLL void SRL_Compare_100x100_ZalesakTest(Scene*);

SRLDLL void SRL_Measure_Scalability(Scene*);