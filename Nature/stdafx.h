// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers
// Windows Header Files:
#include <windows.h>
#include <cstdlib>
#include <iostream>
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glut.h>

// TODO: reference additional headers your program requires here
#include <Nature/Fluid/FluidSimulator/Smoke2D.h>
#include <Nature/Fluid/FluidSimulator/Smoke3D.h>
#include <Nature/Fluid/FluidSimulator/Liquid2D.h>
#include <Nature/Fluid/FluidSimulator/Liquid3D.h>

#include <Principle/Arithmetic/ArrayUtil/ArrayUtil2.h>
#include <Principle/Arithmetic/Arithmetic/CubicEquationSolver_Cardano.h>
#include <Nature/MassSpring/ClothScene/ClothScene.h>

#include <Nature/Fluid/InterfaceTracking/InterfaceIntegrator_SLC.h>