#include <stdafx.h>
#include <Nature/Fluid/CFDsolver/CFDAdvection.h>
#include <Principle/Geometry/MacGrid/MacGridFunctor.h>
#include <Principle/Arithmetic/MathLibrary/Interpolation.h>
#include <Principle/Geometry/Levelset/UpwindDerivative_Functor.h>

using namespace std;
using namespace arithmetic;
using namespace fluidspace;
