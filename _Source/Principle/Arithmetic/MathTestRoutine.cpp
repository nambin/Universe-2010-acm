#include <stdafx.h>
#include <Principle/Arithmetic/MathTestRoutine.h>
#include <Basis/Basis/Array.h>

using namespace std;
using namespace basis;
using namespace arithmetic;

namespace arithmetic {

typedef Vector3<double> Vector3r;
typedef VectorN<double> VectorNr;
typedef Matrix4<double> Matrix4r;
typedef MatrixN<double> MatrixNr;

//typedef float real;
typedef double real;

void CubicEquationSolverTest() {
	CubicEquationSolver_Cardano<real> solver;

	// introduction
	cout << "Cubic Equation Solver" << endl;
	cout << "ax^3 + bx^2 + cx + d = 0" << endl;
	cout << "with a, b, c, d are real, and a is not zero" << endl << endl;
	
	// infinite loop
	while (1) {
		real a,b,c,d;

		// get the co-efficients of x
		cout << "a = ";	cin >> a;
		cout << "b = ";	cin >> b;
		cout << "c = ";	cin >> c;
		cout << "d = ";	cin >> d;
		solver.solve(a, b, c, d);		
		cout << solver;

		for(size_t i=0;i<3;++i) {
			cout << Complex<real>(a) * solver.sol[i] * solver.sol[i] * solver.sol[i] + Complex<real>(b) * solver.sol[i] * solver.sol[i] + Complex<real>(c) * solver.sol[i] + Complex<real>(d) << endl;
		}
		cout << solver.numReal << endl;

		//a * solver.sol[0];
		Complex<real> com;
		a * com;
	}	
}

void SparseMatrixTest() {

	SparseMatrix<double> sparse(3);
	sparse.set_element(0,0,1);

	VectorN<double> v0(3), v1(3);
	mul(v0, sparse, v1);

	SparseColumnLowerFactor<double> fac;
	fac.makePreconditioner_IC0(sparse);

	PCGSolver<double> solver(sparse);

	double res; size_t iteration;
	solver.pcgSolve(v0, v1, res, iteration);

}

void SVDTest() {
	int nDim0 = 6, nDim1 = 7;

	VectorN<real> x(nDim1), b(nDim0);
	MatrixN<real> mat(nDim0, nDim1);	

	for(int i=0;i<nDim0;++i) {
		for(int j=0;j<nDim1;++j) {
			mat(i,j) = rand() / (double)RAND_MAX * 10.0f;
		}		
	}
	for(int j=0;j<nDim0;++j) {
		b(j) = rand() / (double)RAND_MAX * 10.0f;
	}	

	cout << mat << endl;

	SVD<real> svd(mat);

	MatrixN<real> vt(svd.n, svd.n);	vt = transpose(svd.v);
	MatrixN<real> w(svd.n, svd.n);	for(int i=0;i<svd.n;++i) w(i,i) = svd.w(i);
	MatrixN<real> u(svd.m, svd.n);	u = svd.u;

	cout << u << endl << w << endl << vt;

	w = w * vt;
	u = u * w;

	cout << u << endl;

	cout << b << endl;
	svd.solve(x,b);
	b = u * x;
	cout << b << endl;

	MatrixN<real> pseudoInv(1,1);
	svd.getPseudoInverse(pseudoInv);

	x = pseudoInv * b;
	b = u * x;
	cout << b << endl;

//	svd.getRange(mat, 1.e-5);
//	svd.getNullspace(mat, 1.e-5);

//	svd.solve(vec, vec, 1.e-5);
	
}

void VectorMatrixPCATest() {

	int dim0 = 4, dim1 = 10;

	VectorNr vec0, vec1, vec2, vec3;
	MatrixNr mat0, mat1, mat2;
	
	// Vector Initialization	
	vec0.resize(dim0);	vec1.resize(dim0);	vec2.resize(dim0);	vec3.resize(dim0);
	for(int i=0;i<dim0;++i) {
		vec0[i] = rand() / (double)RAND_MAX;
		vec1[i] = rand() / (double)RAND_MAX * 1000.0f;
	}

	vec2 = vec0 + vec1;
	vec3 = vec0 - vec1;

	cout << vec0 << endl;cout << vec1 << endl;cout << vec2 << endl;
	cout << vec3 << endl << endl;

	// Matrix Initialization
	mat0.resize(dim0,dim0);	mat1.resize(dim0,dim1);	mat2.resize(dim1,dim0);
	for(std::size_t i=0;i<mat1.dimRow();++i) {
		for(std::size_t j=0;j<mat1.dimCol();++j) {
			mat1(i,j) = i + j*(i+1) + rand() / (double)RAND_MAX * 10.0f;
		}
	}

	//transpose(mat2, mat1);
	//mul(mat0, mat1, mat2);		// set mat0 = mat1 * mat1^T

	mat2 = transpose(mat1);
	mat0 = mat1 * mat2;

	cout << mat0;
	cout << mat1;
	cout << mat2;

	vec1[0] = 1.0f;	vec1[1] = 2.0f;	vec1[2] = 3.0f;	vec1[3] = 4.0f;

	cout << vec0 << endl;
	cout << vec1 << endl;

	//vec0 = vec0-vec0;

	cout << vec0 << endl;
	cout << vec1 << endl;

	mat0.identity(); mat0 *= 2.0f;

	//vec1 = mat0 * vec0;
	//vec1 = ((double)3.0f) * (mat0 * vec0);
	vec1 = vec1 + (double)3.0f * (mat0 * vec0);
	cout << mat0 << endl;
	cout << vec0 << endl;
	cout << vec1 << endl;
}

};