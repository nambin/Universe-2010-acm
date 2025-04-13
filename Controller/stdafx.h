// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers

// Windows Header Files:
#include <windows.h>

#undef min
#undef max

#include <gl/glew.h>

#include <cstdlib>
#ifdef _WIN32
#include "gl/glut.h"
#else
#include <GLUT/glut.h>
#endif

#include <iostream>
#include <sstream>
#include <string>

// TODO: reference additional headers your program requires here
#include <Basis/Basis/Array2.h>
#include <Basis/Timer/Timer.h>
#include <Basis/DataPath/DataPath.h>

#include <Principle/Arithmetic/MathLibrary/SymmEig.h>
#include <Principle/Arithmetic/Matrix/Matrix4.h>
#include <Principle/Arithmetic/MathTestRoutine.h>
#include <Principle/VirtualSpace/Material/Material.h>

#include <Principle/Geometry/GeometryLibrary/IntersectionTestLibrary.h>
#include <Principle/Geometry/GeometryLibrary/DistanceLibrary.h>
#include <Principle/Geometry/SRL/SRL2D.h>
#include <Principle/Geometry/MeshViewer/MeshViewer.h>
#include <Principle/Geometry/MeshViewer/TriangleMeshViewer.h>
#include <Principle/Geometry/BiTree/QuadTreeViewer.h>
#include <Principle/Geometry/BiTree/OctreeViewer.h>
#include <Principle/Geometry/DistanceFieldLibrary/ADFQuadTree.h>
#include <Principle/Geometry/GeometricPrimitive/Plane.h>
#include <Principle/Geometry/GeometricPrimitive/Sphere.h>
#include <Principle/Geometry/GeometricPrimitive/Circle.h>
#include <Principle/Geometry/bvh/BVH.h>

#include <Illumination/IlluminationHeader/Illumination/Environments.h>
#include <Illumination/IlluminationHeader/Core/ONB.h>
#include <Illumination/IlluminationHeader/Scene/Scene.h>
#include <Illumination/IlluminationHeader/IlluminationTestFunctions.h>
#include <Illumination/RealtimeShader/ShaderViewer/ShaderViewer.h>

#include <Nature/Fluid/Environment/Environments.h>
#include <Nature/Fluid/FluidSimulator/Smoke2D.h>
#include <Nature/Fluid/FluidSimulator/Smoke3D.h>
#include <Nature/Fluid/FluidSimulator/Liquid2D.h>
#include <Nature/Fluid/FluidSimulator/Liquid3D.h>
#include <Nature/Fluid/InterfaceTracking/InterfaceTrackingTest.h>
#include <Nature/Fluid/FluidScene/FluidScene.h>
#include <Nature/Fluid/InterfaceTracking/InterfaceTracking_Viewer.h>
#include <Nature/Fluid/InterfaceTracking/ReinitializationAccuracyTest.h>

#include <Nature/MassSpring/ClothScene/ClothScene.h>
#include <Nature/MassSpring/ClothScene/ClothScene_Objects.h>
#include <Nature/MassSpring/TestFunctions.h>

#include <Nature/RigidBody/Core/StaticRigidBody.h>

#include "Header/TestFunction.h"
#include "Header/TestOpenGLScene.h"

#include "../SRL_Test/SRL_Experiment.h"