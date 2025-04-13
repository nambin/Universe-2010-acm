#pragma once

#include "Scene/Scene.h"
#include "Core/Spectrum.h"
#include "Core/Ray.h"
#include "Core/ONB.h"
#include "Core/Transform.h"
#include "Core/Filter.h"
#include "Core/Camera.h"
#include "Shape/IsoSurface.h"
#include "Shape/IsoSurface_RegularGrid.h"
#include "Shape/Pallelogram.h"
#include "Texture/Texture.h"
#include "Integrator/NormalMapIntegrator.h"
#include "Integrator/DirectLightingIntegrator.h"
#include "Primitive/GeometricPrimitive.h"
#include "Primitive/SimpleAggregate.h"
#include "Primitive/KdTreeAccel.h"
#include "Primitive/GridAccel.h"
#include "Shape/Shape.h"
#include "Shape/Pallelogram.h"
#include "Shape/Sphere.h"
#include "Shape/TriangleMesh.h"
#include "Shape/Cube.h"

namespace illuminationspace {

ILLUMINATIONDLL void test();

ILLUMINATIONDLL void makeScene_CornellBox();
ILLUMINATIONDLL void makeScene_Levelset_ImplicitRendering();
ILLUMINATIONDLL void makeScene_SRL_ImplicitRendering();

ILLUMINATIONDLL void set_CornellBox(std::vector<Primitive> &prims);	
ILLUMINATIONDLL void set_CornellSphere(std::vector<Primitive> &prims);	
ILLUMINATIONDLL void set_Levelset_ImplicitPrimitive(std::vector<Primitive> &prims);	
ILLUMINATIONDLL void set_SRL_ImplicitPrimitive(std::vector<Primitive> &prims);	

ILLUMINATIONDLL void set_BBox(BBox &);
ILLUMINATIONDLL void set_Lights(std::vector<Primitive> &prims, std::vector<Light> &lts);	


};