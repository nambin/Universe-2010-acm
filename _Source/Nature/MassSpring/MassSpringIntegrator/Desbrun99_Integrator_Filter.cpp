#include <stdafx.h>
#include <Nature/MassSpring/MassSpringIntegrator/Desbrun99_Integrator.h>
#include <Principle/Arithmetic/MathLibrary/SVD.h>
#include <Basis/Timer/Timer.h>

using namespace std;
using namespace arithmetic;
using namespace springspace;

void Desbrun99_Integrator::makeFilter(arithmetic::MatrixN<real>& mat, const Array<Spring<real>>& springs, real k) const {
	for(size_t i=0;i<springs.size();++i) {
		const Spring<real>& s = springs[i];

		mat(s.idx0, s.idx1) += k;
		mat(s.idx1, s.idx0) += k;

		mat(s.idx0, s.idx0) -= k;
		mat(s.idx1, s.idx1) -= k;
	}
}

void Desbrun99_Integrator::makeFilter(real dt, const QuadDynamicMesh<real>& mesh) {
	if(mesh.size() == 0)
		basis::throwError("Should initialize QuadMesh first\n");

	basis::Timer timer;

	arithmetic::MatrixN<real> mat;
	mat.init(mesh.size(), mesh.size());

	if(mesh.stretchSpring.size() != 0)	makeFilter(mat, mesh.stretchSpring, mesh.stretchFunc->k);
	if(mesh.shearSpring.size() != 0)	makeFilter(mat, mesh.shearSpring, mesh.shearFunc->k);

	real _mul = dt * dt / mesh.mass;
	for(arithmetic::MatrixN<real>::index_type i=0;i<mat.dimRow();++i) {
		for(arithmetic::MatrixN<real>::index_type j=0;j<mat.dimCol();++j) {
			if(i == j) 	mat(i,i) = 1 - _mul * mat(i,i);
			else 		mat(i,j) = - _mul * mat(i,j);
		}
	}
	
	arithmetic::SVD<real> svd(mat);
	svd.getPseudoInverse(filter);

	timer.printTime("Desbrun99_Integrator::makeFilter() is done");

	if(isFilterValid())	cout << "Filter is Valid" << endl;
	else				cout << "Filter is Not Valid" << endl;
}

bool Desbrun99_Integrator::isFilterValid() const {
	if(!filter.isSymmetricMatrix())
		return false;

	for(MatrixN<real>::index_type i=0;i<filter.dimRow();++i) {
		real sum = 0;
		for(MatrixN<real>::index_type j=0;j<filter.dimCol();++j) 
			sum += filter(i,j);
		if(!IsAlmostSame(sum, real(1)))
			return false;
	}

	return true;
}

void Desbrun99_Integrator::filterForce(Vec3Array<real>& force) const {
	Array<Vector3<real>> forceFiltered(force.size());
	for(Array2<Vector3<real>>::index_type i=0;i<forceFiltered.size();++i) {
		for(Array2<Vector3<real>>::index_type k=0;k<forceFiltered.size();++k) {
			forceFiltered(i) += force(k) * filter(i, k);
		}
	}
	force.set(forceFiltered);
}