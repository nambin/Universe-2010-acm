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

#undef min
#undef max

// TODO: reference additional headers your program requires here

#include <Illumination/IlluminationHeader/Scene/Scene.h>
#include <Illumination/IlluminationHeader/Core/Spectrum.h>
#include <Illumination/IlluminationHeader/Core/Ray.h>
#include <Illumination/IlluminationHeader/Core/ONB.h>
#include <Illumination/IlluminationHeader/Core/Transform.h>
#include <Illumination/IlluminationHeader/Core/Filter.h>
#include <Illumination/IlluminationHeader/Core/Camera.h>
#include <Illumination/IlluminationHeader/Texture/Texture.h>
#include <Illumination/IlluminationHeader/Integrator/NormalMapIntegrator.h>
#include <Illumination/IlluminationHeader/Integrator/DirectLightingIntegrator.h>
#include <Illumination/IlluminationHeader/Integrator/PhotonMapIntegrator.h>
#include <Illumination/IlluminationHeader/Primitive/GeometricPrimitive.h>
#include <Illumination/IlluminationHeader/Primitive/SimpleAggregate.h>
#include <Illumination/IlluminationHeader/Primitive/KdTreeAccel.h>
#include <Illumination/IlluminationHeader/Primitive/GridAccel.h>
#include <Illumination/IlluminationHeader/Shape/Shape.h>
#include <Illumination/IlluminationHeader/Shape/IsoSurface_RegularGrid.h>
#include <Illumination/IlluminationHeader/Shape/Pallelogram.h>
#include <Illumination/IlluminationHeader/Shape/Sphere.h>
#include <Illumination/IlluminationHeader/Shape/TriangleMesh.h>

#include <Illumination/RealtimeShader/ShaderViewer/ShaderViewer.h>
#include <Illumination/RealtimeShader/cgShader/cgShader_VolumeRendering.h>
