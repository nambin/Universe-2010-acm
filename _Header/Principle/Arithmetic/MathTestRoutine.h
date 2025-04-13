#pragma once

#include "../../Basis/Basis/Array.h"
#include "../../Basis/Basis/Array2.h"
#include "../../Basis/Random/RandomGen.h"

#include "ArrayUtil/ArrayUtil2.h"
#include "ArrayUtil/ArrayUtil3.h"
#include "Vector/Vector2.h"
#include "Vector/Vector3.h"
#include "Vector/Vector4.h"
#include "Vector/VectorN.h"
#include "Matrix/Matrix3.h"
#include "Matrix/Matrix4.h"
#include "Matrix/MatrixN.h"
#include "Matrix/FixedSparseMatrix.h"
#include "Matrix/PCGSolver.h"
#include "MathLibrary/SymmEig.h"
#include "MathLibrary/PCA.h"
#include "MathLibrary/SVD.h"
#include "Arithmetic/CubicEquationSolver_Cardano.h"
#include "../Principle.h"

namespace arithmetic {

void PRINCIPLEDLL CubicEquationSolverTest();
void PRINCIPLEDLL VectorMatrixPCATest();
void PRINCIPLEDLL SVDTest();
void PRINCIPLEDLL SparseMatrixTest();

};		// namespace arithmetic