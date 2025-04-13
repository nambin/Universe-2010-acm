// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#define WIN32_LEAN_AND_MEAN             // Exclude rarely-used stuff from Windows headers
// Windows Header Files:
#include <windows.h>

// TODO: reference additional headers your program requires here
#include <Principle/Arithmetic/L_BFGS_B/ap.h>
#include <Principle/Arithmetic/L_BFGS_B/L_BFGS_B.h>
#include <Principle/Arithmetic/MathLibrary/PCA.h>

#include <Basis/Color/RGBA.h>
#include <Principle/VirtualSpace/Light/Light.h>
#include <Principle/VirtualSpace/Material/Material.h>
#include <Principle/VirtualSpace/OpenGLState/OpenGLState.h>
#include <Principle/VirtualSpace/Property/MatLight.h>
#include <Principle/VirtualSpace/Scene/Scene.h>
#include <Principle/VirtualSpace/DrawLibrary/glDrawLibrary.h>

#include <Principle/Geometry/BiTree/QuadTreeViewer.h>
#include <Principle/Geometry/MeshViewer/MeshViewer.h>